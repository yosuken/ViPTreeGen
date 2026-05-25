#
#  duckdb_util.rb / ViPTreeGen - helper for DuckDB CLI access
#
#    Copyright: 2026 (C) Yosuke Nishimura
#    License: MIT license
#
#  Provides a thin wrapper around the `duckdb` CLI so that the pipeline can
#  consolidate intermediate data and logs into a single DuckDB database per run.
#
#  Usage:
#    ENV["DUCKDB_PATH"] = "#{Odir}/run.duckdb"
#    DuckDBUtil.exec_sql(DuckDBUtil::SCHEMA_SQL)
#    DuckDBUtil.copy_in("sequences", "tmp.tsv", columns: %w[seq_id kind length ord])
#    DuckDBUtil.query_each("SELECT seq_id, length FROM sequences ORDER BY ord") { |line| ... }
#

module DuckDBUtil
  module_function

  SCHEMA_SQL = <<~SQL
    CREATE TABLE IF NOT EXISTS run_metadata (
      key   VARCHAR PRIMARY KEY,
      value VARCHAR
    );

    -- seq: the nucleotide sequence itself, stored so that a finished run.duckdb
    -- is sufficient to seed a downstream with-reference-mode run (--ref-duckdb)
    -- without needing the original FASTA file alongside it. New in schema 2.0.
    CREATE TABLE IF NOT EXISTS sequences (
      seq_id VARCHAR NOT NULL,
      kind   VARCHAR NOT NULL,
      length INTEGER NOT NULL,
      ord    INTEGER NOT NULL,
      seq    VARCHAR NOT NULL,
      PRIMARY KEY (seq_id, kind)
    );

    -- kind: 'node' = a query genome (db = all.fasta); 'input' = self-blast (db = single seq), 2D mode only
    CREATE TABLE IF NOT EXISTS summary_pre (
      kind               VARCHAR NOT NULL,
      node               VARCHAR,
      que_score_rank     INTEGER,
      sub_score_rank     INTEGER,
      que                VARCHAR,
      sub                VARCHAR,
      que_score          INTEGER,
      sub_score          INTEGER,
      que_len            INTEGER,
      sub_len            INTEGER,
      que_len_in_hit     INTEGER,
      sub_len_in_hit     INTEGER,
      mean_idt_of_que    DOUBLE,
      mean_idt_of_sub    DOUBLE,
      pct_que_len_in_hit DOUBLE,
      pct_sub_len_in_hit DOUBLE
    );

    CREATE TABLE IF NOT EXISTS self_scores (
      seq_id     VARCHAR PRIMARY KEY,
      self_score DOUBLE
    );

    CREATE TABLE IF NOT EXISTS summary_tsv (
      node     VARCHAR,
      target   VARCHAR,
      length   INTEGER,
      score    DOUBLE,
      sg       DOUBLE,
      mean_idt DOUBLE,
      pct_len  DOUBLE
    );

    CREATE TABLE IF NOT EXISTS logs (
      path      VARCHAR,
      source    VARCHAR,
      seq_id    VARCHAR,
      split_idx INTEGER,
      ts        TIMESTAMP,
      content   VARCHAR
    );

    -- resolved path + reported version of each external tool actually used by
    -- this run (engine binaries, duckdb, ruby, R, parallel, the bundled Rust
    -- binary). Populated by InitDuckDB from the CLI at run start.
    CREATE TABLE IF NOT EXISTS tools (
      name    VARCHAR PRIMARY KEY,
      path    VARCHAR,
      version VARCHAR
    );
  SQL

  def db_path
    ENV["DUCKDB_PATH"] or raise("DUCKDB_PATH env var not set")
  end

  def duckdb_bin
    ENV["DUCKDB_BIN"] || "duckdb"
  end

  # Execute one or more SQL statements via stdin. Returns stdout as a String.
  # Raises if the duckdb process exits non-zero.
  def exec_sql(sql, db: db_path)
    out = IO.popen([duckdb_bin, db], "r+") do |io|
      io.write(sql)
      io.close_write
      io.read
    end
    unless $?.success?
      raise "DuckDB SQL failed (status=#{$?.exitstatus}):\n--- SQL ---\n#{sql}\n--- stdout ---\n#{out}"
    end
    out
  end

  # Stream the result rows of a query as tab-separated lines (no header).
  # Opens the database read-only so concurrent parallel workers (e.g. via GNU parallel
  # in 02-3) can share access -- DuckDB enforces a single-writer lock otherwise.
  def query_each(sql, db: db_path)
    IO.popen([duckdb_bin, "-readonly", "-list", "-noheader", "-separator", "\t", db, sql], "r") do |io|
      io.each_line { |line| yield line.chomp }
    end
    raise "DuckDB query failed: #{sql}" unless $?.success?
  end

  # Bulk-load a tab-separated file into a table.
  def copy_in(table, tsv_path, columns: nil, db: db_path)
    cols = columns ? " (#{columns.join(',')})" : ""
    sql  = "COPY #{table}#{cols} FROM '#{tsv_path}' (FORMAT CSV, DELIMITER E'\\t', HEADER false);"
    exec_sql(sql, db: db)
  end

  # Ingest a log file's entire content as a single row in `logs`.
  # `read_text` is a table function returning columns (filename, content, size, last_modified);
  # we SELECT its `content` column as the log body.
  def ingest_log(path, source:, seq_id: nil, split_idx: nil, db: db_path)
    return unless File.exist?(path)
    seq_lit   = seq_id    ? sql_str(seq_id) : "NULL"
    split_lit = split_idx ? split_idx.to_i.to_s : "NULL"
    sql = "INSERT INTO logs SELECT #{sql_str(path)}, #{sql_str(source)}, " \
          "#{seq_lit}, #{split_lit}, now(), content FROM read_text(#{sql_str(path)});"
    exec_sql(sql, db: db)
  end

  # Escape a string for SQL literal (single-quoted, doubled internal quotes).
  def sql_str(s)
    "'#{s.to_s.gsub("'", "''")}'"
  end
end
