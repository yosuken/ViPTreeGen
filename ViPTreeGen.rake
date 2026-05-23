#
#  ViPTreeGen.rake - a command line tool for viral proteomic tree generation
#
#    Copyright: 2017 (C) Yosuke Nishimura (yosuke@kuicr.kyoto-u.ac.jp)
#    License: MIT license
#    Initial version: 2017-02-07
#

require "time"
require_relative "script/duckdb_util"


# {{{ procedures
WriteBatch  = lambda do |outs, jdir, t|
	outs.each_slice(10000).with_index(1){ |ls, idx|
		open("#{jdir}/#{t.name.split(".")[1..-1]*"."}.#{idx}", "w"){ |fjob|
			fjob.puts ls
		}
	}
end

RunBatch    = lambda do |jdir, queue, nthreads, mem, wtime, ncpus|
	# [TODO] queue validation
	Dir["#{jdir}/*"].sort_by{ |fin| fin.split(".")[-1].to_i }.each{ |fin|
		if queue != ""
			raise("`--queue #{queue}': invalid queue") unless %w|JP1 JP4 JP10 cdb|.include?(queue)
			sh "qsubarraywww -q #{queue} -l ncpus=#{nthreads} -l mem=#{mem}gb -l walltime=#{wtime} #{fin}"
		elsif ncpus != ""
			raise("`--ncpus #{ncpus}': not an integer") if ncpus !~ /^\d+$/
			sh "parallel --jobs #{ncpus} <#{fin}"
		else
			sh "sh #{fin}"
		end
	}
end

## --resume support: each Rake task records `step_done:<task.name>` in run_metadata when
## it completes; on a subsequent run with --resume those tasks are skipped. Detection is
## a single SELECT; marking is idempotent (INSERT OR REPLACE on the metadata PK).
StepDone = lambda do |task_name|
	## Cannot use DuckDBUtil.query_each before InitDuckDB has created the DB on a fresh
	## run, so guard on file existence.
	db = ENV["DUCKDB_PATH"]
	return false unless db && File.file?(db)
	done = false
	DuckDBUtil.query_each(
		"SELECT 1 FROM run_metadata WHERE key='step_done:#{task_name}' LIMIT 1"
	) { |_| done = true }
	done
end

MarkStepDone = lambda do |task_name|
	DuckDBUtil.exec_sql(
		"INSERT OR REPLACE INTO run_metadata VALUES " \
		"('step_done:#{task_name}', #{DuckDBUtil.sql_str(Time.now.iso8601)});"
	)
end

PrintStatus = lambda do |current, total, status, t|
	puts ""
	puts "\e[1;32m===== #{Time.now}\e[0m"
	puts "\e[1;32m===== step #{current} / #{total} (#{t.name}) -- #{status}\e[0m"
	puts ""
	$stdout.flush
end

InitDuckDB = lambda do |mode:|
	## (re)create schema and seed run_metadata. Called from 01-1 (both normal and 2D).
	DuckDBUtil.exec_sql(DuckDBUtil::SCHEMA_SQL)
	rows = [
		["version",        "2.0.0-dev"],
		["schema_version", "2.0"],         ## bumped when `sequences.seq` was added; required by --ref-duckdb
		["mode",           mode],          ## "normal", "2D", or "ref"
		["search_mode",    Mode],          ## "tblastx" or "mmseqs-tblastx"
		["ref_mode",       RefMode.to_s],
		["ref_duckdb_path", RefMode ? RefDuckDB : ""],
		["cutlen",      Cutlen.to_s],
		["dbsize",      DBsize.to_s],
		["matrix",      Matrix.to_s],
		["evalue",      Evalue.to_s],
		["tree_method", TreeMethod.to_s],
		["start_time",  Time.now.iso8601],
	]
	sql = "DELETE FROM run_metadata;\n" + rows.map { |k, v|
		"INSERT INTO run_metadata VALUES (#{DuckDBUtil.sql_str(k)}, #{DuckDBUtil.sql_str(v)});"
	}.join("\n")
	DuckDBUtil.exec_sql(sql)
end

FinalLogCleanup = lambda do
	## Ingest per-node tblastx.log files into the `logs` table, then delete them.
	## Also delete the makeblastdb.log (already ingested by 01-1).
	## Top-level log/*.log (user-visible pipeline log) is left in place.
	tblastx_logs = Dir["#{Odir}/node/*/blast/tblastx.log"] + Dir["#{Odir}/input/*/blast/tblastx.log"]
	unless tblastx_logs.empty?
		sql = +"BEGIN;\n"
		tblastx_logs.each do |path|
			node = path.split("/")[-3]
			kind = path.include?("/input/") ? "input" : "node"
			seq_lit = DuckDBUtil.sql_str(node)
			path_lit = DuckDBUtil.sql_str(path)
			sql << "INSERT INTO logs SELECT #{path_lit}, 'tblastx', #{seq_lit}, NULL, now(), content FROM read_text(#{path_lit});\n"
		end
		sql << "COMMIT;\n"
		DuckDBUtil.exec_sql(sql)
	end
	sh "find #{Odir}/node #{Odir}/input -name 'tblastx.log' -delete 2>/dev/null || true"
	sh "rm -f #{Odir}/cat/all/all.fasta.makeblastdb.log"

	## mmseqs per-stage logs ingest into `logs` (source = "mmseqs:<stage>"), then delete.
	## Files come from 01-1 (createdb logs) and 01-2 (search/convertalis logs).
	if UsesMmseqs
		mmseqs_logs = Dir["#{Odir}/cat/all/mmseqs_*.log"] +
		              Dir["#{Odir}/cat/all/*.fasta.createdb.log"]
		unless mmseqs_logs.empty?
			sql = +"BEGIN;\n"
			mmseqs_logs.each do |path|
				base  = File.basename(path)
				## strip "mmseqs_" prefix and ".log" suffix; drop ".fasta" middle segment
				stage = base.sub(/\.log\z/, "").sub(/\Ammseqs_/, "").sub(/\.fasta\./, ".")
				source = "mmseqs:#{stage}"  # e.g. "mmseqs:search.node", "mmseqs:input.createdb"
				path_lit   = DuckDBUtil.sql_str(path)
				source_lit = DuckDBUtil.sql_str(source)
				sql << "INSERT INTO logs SELECT #{path_lit}, #{source_lit}, NULL, NULL, now(), content FROM read_text(#{path_lit});\n"
			end
			sql << "COMMIT;\n"
			DuckDBUtil.exec_sql(sql)
		end
		## mmseqs artifacts (DB / index / tmp / logs)
		sh "rm -rf #{Odir}/cat/all/all.fasta.mmdb* #{Odir}/cat/all/all.fasta.query.mmdb* #{Odir}/cat/all/all.fasta.input.mmdb* #{Odir}/mmseqs_tmp 2>/dev/null || true"
		sh "rm -f #{Odir}/cat/all/mmseqs_*.log #{Odir}/cat/all/*.fasta.createdb.log 2>/dev/null || true"
		sh "find #{Odir}/input -name '*.mmdb*' -delete 2>/dev/null || true"
	end
end

## Parse a FASTA string, normalize sequence labels, and validate per-entry
## constraints. Returns an Array of [normalized_label, sequence] pairs in
## FASTA order. Raises on format errors / too-short sequences / within-file
## duplicate labels. Prints "### [!] X is changed to Y" for each renamed entry.
##
## Label normalization rules (applied in order):
##   - take the first whitespace-delimited token of the header line
##   - replace any char outside [a-zA-Z0-9.\-_] with `_`
##   - strip leading / trailing `.` `-` `_`
##   - collapse consecutive `.` `-` `_` to a single char
##
## Reasons: tblastx/tree generation (ape/phangorn bionj) is sensitive to
## leading/trailing dots and runs of separators.
ParseFastaEntries = lambda do |fasta_str, fasname:, min_seqs: nil|
	raise("\e[1;31mError:\e[0m #{fasname} might not be in FASTA format.") if fasta_str[0] != ">"
	ents = fasta_str.split(/^>/)[1..-1] || []
	raise("\e[1;31mError:\e[0m #{fasname} should include at least #{min_seqs} seqeuences.") if min_seqs && ents.size < min_seqs

	not_allowed_pattern     = /[^a-zA-Z0-9\-\.\_]/
	not_allowed_pattern_pre = /^[\.\-\_]+/
	not_allowed_pattern_suf = /[\.\-\_]+$/
	not_allowed_pattern_dot = /(\.)[\.]+/
	not_allowed_pattern_hyp = /(\-)[\-]+/
	not_allowed_pattern_und = /(\_)[\_]+/

	seen   = {}
	result = []
	ents.each do |ent|
		lab, *seq = ent.split("\n")
		seq  = seq.join.gsub(/\s/, "")
		len  = seq.size

		_lab = lab.split(/\s+/)[0]
		lab  = _lab.gsub(not_allowed_pattern, "_")
		lab  = lab.gsub(not_allowed_pattern_pre, "")
		          .gsub(not_allowed_pattern_suf, "")
		          .gsub(not_allowed_pattern_dot, '\1')
		          .gsub(not_allowed_pattern_hyp, '\1')
		          .gsub(not_allowed_pattern_und, '\1')

		puts "### [!] #{_lab} is changed to #{lab}" if _lab != lab

		raise("\e[1;31mError:\e[0m '#{lab}' is too short (length < 100 nt) included in #{fasname}.") if len < 100
		raise("\e[1;31mError:\e[0m sequence name is not uniq. '#{lab}' is found multiple times in #{fasname}.") if seen[lab]

		seen[lab] = true
		result << [lab, seq]
	end
	result
end

IngestSequences = lambda do |labs_in_order, lab2seq, kind:, ord_offset: 0|
	## bulk-insert into sequences table via tmp TSV. `seq` is the nucleotide
	## string itself (no whitespace — ParseFastaEntries stripped it), so a 4-col
	## TSV row stays unambiguous and tab-safe. `ord_offset` lets ref-mode append
	## input sequences after pre-imported ref entries without ord collisions.
	require "fileutils"
	tdir = "#{Odir}/tmp"; FileUtils.mkdir_p(tdir)
	tsv  = "#{tdir}/sequences.#{kind}.tsv"
	open(tsv, "w") { |f|
		labs_in_order.each_with_index { |lab, idx|
			seq = lab2seq[lab]
			next unless seq
			f.puts [lab, kind, seq.size, idx + ord_offset, seq].join("\t")
		}
	}
	DuckDBUtil.copy_in("sequences", tsv, columns: %w[seq_id kind length ord seq])
	File.delete(tsv)
end

## Build the command that creates a search-engine database from a FASTA file.
## In `tblastx` mode this is `makeblastdb` against the FASTA in place; in `mmseqs`
## mode it is `mmseqs createdb` (optionally followed by `mmseqs createindex` for
## the large target DB). The lambda returns one shell command string suitable
## for `sh`. `db_path` is the prefix that subsequent search commands will pass
## as `-db` / target: in tblastx mode this equals `fasta`, in mmseqs mode it has
## a `.mmdb` suffix so it does not collide with the FASTA file.
## Build a BLAST+ search command for one (query, db, out, log) tuple. Used by 01-1 /
## 01-1.2D / 01-1.ref.prep to populate per-node batch files. mmseqs-* modes do their
## search in 01-2 (single global call) and never enter this helper.
SearchCmd = lambda do |query:, db:, out:, log:|
	case Mode
	when "tblastx"
		"tblastx -dbsize #{DBsize} -matrix #{Matrix} -max_target_seqs #{Max_target_seqs} " \
		"-num_threads #{Nthreads} -evalue #{Evalue} -outfmt 6 -db #{db} -query #{query} -out #{out} 2>#{log}"
	when "blastn"
		"blastn -dbsize #{DBsize} -max_target_seqs #{Max_target_seqs} " \
		"-num_threads #{Nthreads} -evalue #{Evalue} -outfmt 6 -db #{db} -query #{query} -out #{out} 2>#{log}"
	else
		raise "SearchCmd: unsupported Mode '#{Mode}' for per-batch search (mmseqs-* modes use one global call in 01-2)"
	end
end

BuildDBCmd = lambda do |mode:, fasta:, db_path:, log:, create_index: false|
	## All search modes here run against a nucleotide DB (tblastx does on-the-fly
	## translation; mmseqs-tblastx similarly. The DB itself stays nucleotide.)
	case mode
	when "tblastx", "blastn"
		"makeblastdb -dbtype nucl -in #{fasta} -out #{db_path} -title #{File.basename(fasta)} 2>#{log}"
	when "mmseqs-tblastx"
		idx_tmp = "#{File.dirname(db_path)}/mmseqs_idx_tmp.#{File.basename(db_path)}"
		cmds = ["mmseqs createdb #{fasta} #{db_path} >#{log} 2>&1"]
		if create_index
			cmds << "mkdir -p #{idx_tmp}"
			cmds << "mmseqs createindex #{db_path} #{idx_tmp} --search-type 4 >>#{log} 2>&1"
			cmds << "rm -rf #{idx_tmp}"
		end
		cmds.join(" && ")
	else
		raise "BuildDBCmd: unknown mode '#{mode}'"
	end
end

RustBin = lambda do
	## Locate the bundled Rust binary `viptreegen-summary-pre`.
	##   1. $VIPTREEGEN_SUMMARY_PRE override
	##   2. `which viptreegen-summary-pre`           (Bioconda / system install)
	##   3. ./rust/target/release/viptreegen-summary-pre  (dev mode after `cargo build --release`)
	if env = ENV["VIPTREEGEN_SUMMARY_PRE"]
		return env
	end
	path_bin = `which viptreegen-summary-pre 2>/dev/null`.chomp
	return path_bin unless path_bin.empty?
	dev_bin = File.expand_path("rust/target/release/viptreegen-summary-pre", __dir__)
	return dev_bin if File.executable?(dev_bin)
	raise "\e[1;31mError:\e[0m viptreegen-summary-pre not found. " \
	      "Install via Bioconda, or run: cargo build --release --manifest-path #{File.expand_path('rust/Cargo.toml', __dir__)}"
end

CheckVersion = lambda do |commands|
	commands.each{ |command|
		str = case command
					when "ruby"
						%|ruby --version 2>&1|
					when "makeblastdb"
						%|makeblastdb -version 2>&1|
					when "tblastx"
						%|tblastx -version 2>&1|
					when "blastn"
						%|blastn -version 2>&1|
					when "R"
						%|LANG=C R --version 2>&1|
					when "ape"
						%|LANG=C R --quiet --no-save --no-restore -e "packageVersion('ape')" 2>&1|
					when "phangorn"
						%|LANG=C R --quiet --no-save --no-restore -e "packageVersion('phangorn')" 2>&1|
					when "parallel"
						%{LANG=C parallel --version 2>&1 |head -n 1}
					when "duckdb"
						%|duckdb --version 2>&1|
					when "mmseqs"
						%|mmseqs version 2>&1|
					when "viptreegen-summary-pre"
						%{#{RustBin.call} --help 2>&1 | head -n 2}
					end
		puts ""
		puts "\e[1;32m===== check version: #{command}\e[0m"
		puts ""
		puts "$ #{str}"
		### run
		puts `#{str}`
		### flush
		$stdout.flush
	}
end
# }}} procedures


# {{{ default (run all tasks)
task :default do
	### define shared tasks
	tasks = %w|01-2.tblastx 02-3.make_summary_pre 02-4.make_self_tblastx|

	### add specific tasks
	if ENV["twoD"] != ""
		## 2D mode (mutex with ref mode; CLI already enforced)
		tasks = %w|01-1.2D.prep_for_tblastx| + tasks + %w|02-5.2D.make_summary_tsv 03-1.2D.make_matrix|
	elsif ENV["ref_duckdb"].to_s != ""
		## with-reference mode: 01-1.ref.prep replaces 01-1.prep_for_tblastx
		## (final tree/matrix tasks are identical to normal mode — the matrix is over input ∪ ref)
		if ENV["notree"] == ""
			tasks = %w|01-1.ref.prep| + tasks + %w|02-5.make_summary_tsv 03-1.make_matrix 03-2.matrix_to_nj|
		else
			tasks = %w|01-1.ref.prep| + tasks + %w|02-5.make_summary_tsv 03-1.make_matrix|
		end
	else
		## normal mode
		if ENV["notree"] == "" ## make tree
			tasks = %w|01-1.prep_for_tblastx| + tasks + %w|02-5.make_summary_tsv 03-1.make_matrix 03-2.matrix_to_nj|
		else ## only generate matrix
			tasks = %w|01-1.prep_for_tblastx| + tasks + %w|02-5.make_summary_tsv 03-1.make_matrix|
		end
	end

	### constants
	Odir     = ENV["dir"]
	Fin      = ENV["fin"]                  ## input file when 2D & normal mode
	Fin_q    = ENV["twoD"]                 ## query file only when 2D mode
	Flen     = "#{Odir}/cat/all/all.len"   ## create when 2D & normal mode
	Flen_q   = "#{Odir}/cat/all/query.len" ## create only when 2D mode
	Flen_i   = "#{Odir}/cat/all/input.len" ## create only when 2D mode
	Fa       = "#{Odir}/cat/all/all.fasta" ## create when 2D & normal mode
	Ddb      = "#{Odir}/run.duckdb"        ## single DuckDB file collecting tmp data and logs
	ENV["DUCKDB_PATH"] = Ddb               ## expose to script/*.rb

	Cutlen   = ENV["cutlen"].to_i ## default: 100,000
	DBsize   = ENV["dbsize"]      ## default: 200,000,000
	Matrix   = ENV["matrix"]      ## default: BLOSUM45
	Evalue   = ENV["evalue"]      ## default: 1e-2
	## Search engine + algorithm pairing.
	##   - tblastx        : NCBI BLAST+ tblastx        (translated 6-frame, protein-level scoring) -- proteomic tree
	##   - mmseqs-tblastx : MMseqs2 --search-type 4    (translated 6-frame, protein-level scoring) -- proteomic tree (default)
	##   - blastn         : NCBI BLAST+ blastn         (nucleotide-vs-nucleotide)                  -- nucleotide tree (DiGAlign backend)
	## ('mmseqs-blastn' was evaluated and removed; see doc/mmseqs-blastn.md.)
	Mode     = ENV["mode"] || "mmseqs-tblastx"
	valid_modes = %w|tblastx mmseqs-tblastx blastn|
	raise("`--mode #{Mode}': must be one of #{valid_modes.join(', ')}") unless valid_modes.include?(Mode)
	## Convenience flag: true for mmseqs-tblastx (single global search in 01-2 instead of per-(node,split) batches).
	UsesMmseqs = (Mode == "mmseqs-tblastx")
	## with-reference mode: previous-run `run.duckdb` providing ref sequences, self_scores, summary_tsv.
	RefDuckDB = ENV["ref_duckdb"].to_s
	RefMode   = !RefDuckDB.empty?
	## --resume: continue from the last completed step (step_done markers in run_metadata).
	Resume    = ENV["resume"].to_s == "true"
	## In mmseqs modes the search DB uses a different prefix; in BLAST+ modes it equals the FASTA path.
	DB       = UsesMmseqs ? "#{Fa}.mmdb" : Fa
	## mmseqs target index memory cap (passed to mmseqs `--split-memory-limit`).
	## Default 12G; user can override via --mmseqs-split-memory-limit.
	MmseqsSplitMem = ENV["mmseqs_split_memory_limit"].to_s.empty? ? "12G" : ENV["mmseqs_split_memory_limit"]
	Nthreads = "1".to_i
	Mem      = Nthreads * 12
	Qname    = ENV["queue"]||""
	Wtime    = ENV["wtime"]||"24:00:00"
	Ncpus    = ENV["ncpus"]||""

	Max_target_seqs = 1_000_000

	TreeMethod = ENV["method"]||"bionj"

	### check version
	engine_cmds = if UsesMmseqs
		%w|mmseqs|
	elsif Mode == "tblastx"
		%w|tblastx makeblastdb|
	else  # blastn
		%w|blastn makeblastdb|
	end
	commands    = engine_cmds + %w|duckdb viptreegen-summary-pre R ape phangorn ruby|
	commands += %w|parallel| if Ncpus != ""
	CheckVersion.call(commands)

	### resume preflight: verify previous-run parameters match, and clean up partial state
	### from whatever step was killed mid-flight. Skipped on a fresh run.
	if Resume
		raise "\e[1;31mError:\e[0m --resume specified but #{Ddb} not found" unless File.file?(Ddb)
		prev = {}
		DuckDBUtil.query_each("SELECT key, value FROM run_metadata") do |line|
			k, v = line.split("\t", 2)
			prev[k] = v
		end

		## parameter integrity: compare critical knobs against the previous run.
		twoD_now = (ENV["twoD"] || "") != ""
		run_mode = if twoD_now then "2D"
		           elsif RefMode then "ref"
		           else "normal"
		           end
		checks = {
			"schema_version"   => "2.0",
			"mode"             => run_mode,
			"search_mode"      => Mode,
			"cutlen"           => Cutlen.to_s,
			"dbsize"           => DBsize.to_s,
			"matrix"           => Matrix.to_s,
			"evalue"           => Evalue.to_s,
			"ref_mode"         => RefMode.to_s,
			"ref_duckdb_path"  => RefMode ? RefDuckDB : "",
		}
		diffs = checks.reject { |k, v| prev[k] == v || (prev[k].nil? && v.to_s.empty?) }
		unless diffs.empty?
			msg = diffs.map { |k, v| "  #{k}: prev=#{prev[k].inspect}, current=#{v.inspect}" }.join("\n")
			raise "\e[1;31mError:\e[0m --resume rejected: parameters differ from the prior run:\n#{msg}\n" \
			      "Either re-invoke with the original parameters, or start a fresh run in a new output dir."
		end
		puts "\e[1;36m### resume preflight OK -- parameters match the prior run\e[0m"

		## cleanup partial mid-step state. Steps that completed atomically (in a SQL
		## transaction) are unaffected; this handles the on-disk artifacts of the
		## search step (01-2) and the mmseqs tmpdir.
		sh "rm -rf #{Odir}/mmseqs_tmp 2>/dev/null || true"
		unless StepDone.call("01-2.tblastx")
			## tblastx mode: partial GNU parallel batches may have written some tblastx.out
			## files. Re-run 01-2 produces them fresh; pre-delete to avoid mixing partial
			## results from the killed run with new output.
			sh "find #{Odir}/node #{Odir}/input -type f -name 'tblastx.out' -delete 2>/dev/null || true"
		end
	end

	### run
	NumStep  = tasks.size
	tasks.each.with_index(1){ |task, idx|
		Rake::Task[task].invoke(idx)
	}
end
# }}} default (run all tasks)


# {{{ tasks 01
desc "01-1.2D.prep_for_tblastx"
task "01-1.2D.prep_for_tblastx", ["step"] do |t, args|
	if Resume && StepDone.call(t.name)
		PrintStatus.call(args.step, NumStep, "SKIP (resume)", t)
		next
	end
	PrintStatus.call(args.step, NumStep, "START", t)
	jdir     = "#{Odir}/batch/#{t.name.split(".")[0]}"; mkdir_p jdir
	odir     = "#{Odir}/cat/all"; mkdir_p odir
	log      = "#{odir}/all.fasta.makeblastdb.log"
	outs     = []
	puts ""

	## initialize DuckDB and record run metadata
	InitDuckDB.call(mode: "2D")

	## validate and make copy of input fasta
	fasta_str_q = IO.read(Fin_q)
	fasta_str_i = IO.read(Fin)

	entries_q = ParseFastaEntries.call(fasta_str_q, fasname: "query FASTA file")
	entries_i = ParseFastaEntries.call(fasta_str_i, fasname: "input FASTA file", min_seqs: 3)

	lab2seq_q = {}
	lab2seq_i = {}
	uniqlab   = {}
	[[entries_q, lab2seq_q], [entries_i, lab2seq_i]].each{ |entries, lab2seq|
		entries.each{ |lab, seq|
			lab2seq[lab] = seq
			uniqlab[lab] = 1
		}
	}
	### write all.fasta & all.len
	open(Fa, "w"){ |fall|
		open(Flen, "w"){ |flen_a|
			uniqlab.each_key{ |lab|
				if q = lab2seq_q[lab] and i = lab2seq_i[lab]
					### detect same name but different sequence entry
					raise("\e[1;31mError:\e[0m same name entries are included in both input FASTA and query FASTA, but the sequences differ.") if q != i
				end
				seq = lab2seq_q[lab]||lab2seq_i[lab]
				fall.puts [">"+lab, seq.scan(/.{1,70}/)]
				flen_a.puts [lab, seq.size]*"\t"
			}
		}
	}
	### write input.len & query.len
	open(Flen_i, "w"){ |flen_i|
		open(Flen_q, "w"){ |flen_q|
			[lab2seq_q, lab2seq_i].zip([flen_q, flen_i]){ |lab2seq, flen|
				lab2seq.each{ |lab, seq|
					flen.puts [lab, seq.size]*"\t"
				}
			}
		}
	}

	### write node/input fasta
	[lab2seq_q, lab2seq_i].zip(%w|node input|){ |lab2seq, type|
		lab2seq.each{ |lab, seq|
			len   = seq.size
			n0dir = "#{Odir}/#{type}/#{lab}/seq";   mkdir_p n0dir
			n1dir = "#{Odir}/#{type}/#{lab}/blast"; mkdir_p n1dir

			fque   = "#{n0dir}/#{lab}.fasta"
			self_db = (UsesMmseqs) ? "#{fque}.mmdb" : fque
			open(fque, "w"){ |_fque|
				_fque.puts [">"+lab, seq.scan(/.{1,70}/)]
			}

			if !UsesMmseqs && type == "input"
				## per-input self-DB; needed in BLAST+ modes (tblastx / blastn) for input self-blast.
				## mmseqs uses one bulk input-vs-input search instead.
				sh BuildDBCmd.call(mode: Mode, fasta: fque, db_path: self_db,
				                   log: "#{fque}.makeblastdb.log", create_index: false)
			end

			## mmseqs runs one global search in 01-2 -- skip per-(node, split) batch generation.
			next if UsesMmseqs

			if len > Cutlen
				n0dir_s = "#{n0dir}/split"; mkdir_p n0dir_s

				idx = 0
				while seq and seq.size > 0
					idx += 1
					subseq, seq = seq[0, Cutlen], seq[Cutlen..-1]
					fspt = "#{n0dir_s}/#{idx}.fasta"
					open(fspt, "w"){ |_fspt|
						_fspt.puts [">#{idx}", subseq.scan(/.{1,70}/)]
					}

					n1dir_s = "#{n1dir}/split/#{idx}"; mkdir_p n1dir_s
					_out  = "#{n1dir_s}/tblastx.out"
					_log  = "#{n1dir_s}/tblastx.log"
					_db   = (type == "input") ? self_db : DB
					outs << SearchCmd.call(query: fspt, db: _db, out: _out, log: _log)
				end
			else
				_out = "#{n1dir}/tblastx.out"
				_log = "#{n1dir}/tblastx.log"
				_db  = (type == "input") ? self_db : DB
				outs << SearchCmd.call(query: fque, db: _db, out: _out, log: _log)
			end
		}
	}

	## write batch file
	WriteBatch.call(outs, jdir, t)

	## build target DB for all.fasta (createindex for mmseqs to amortize per-query lookups)
	sh BuildDBCmd.call(mode: Mode, fasta: Fa, db_path: DB, log: log, create_index: (UsesMmseqs))

	## mmseqs 2D mode: build query.fasta / input.fasta DBs for the two separate searches.
	## (tblastx mode uses per-input self-DBs created inline above; mmseqs path replaces
	## that pattern with one input-vs-input bulk search.)
	if UsesMmseqs
		query_fa = "#{Odir}/cat/all/query.fasta"
		input_fa = "#{Odir}/cat/all/input.fasta"
		open(query_fa, "w") { |f| lab2seq_q.each { |l, s| f.puts ">"+l, s.scan(/.{1,70}/) } }
		open(input_fa, "w") { |f| lab2seq_i.each { |l, s| f.puts ">"+l, s.scan(/.{1,70}/) } }
		sh BuildDBCmd.call(mode: Mode, fasta: query_fa, db_path: "#{Fa}.query.mmdb",
		                   log: "#{query_fa}.createdb.log", create_index: false)
		sh BuildDBCmd.call(mode: Mode, fasta: input_fa, db_path: "#{Fa}.input.mmdb",
		                   log: "#{input_fa}.createdb.log", create_index: true)
	end

	## populate DuckDB: sequences (query = kind 'node', input = kind 'input') + makeblastdb.log
	IngestSequences.call(lab2seq_q.keys, lab2seq_q, kind: "node")
	IngestSequences.call(lab2seq_i.keys, lab2seq_i, kind: "input")
	DuckDBUtil.ingest_log(log, source: (UsesMmseqs ? "mmseqs:all.createdb" : "makeblastdb"))
	MarkStepDone.call(t.name)
end
desc "01-1.ref.prep"
task "01-1.ref.prep", ["step"] do |t, args|
	## with-reference-mode prep:
	##   - read ref-duckdb's run_metadata to verify schema_version >= 2.0 and mode == 'normal'
	##   - ParseFastaEntries on input.fasta, detect ID collisions with ref
	##   - InitDuckDB (mode='ref') -> new run.duckdb
	##   - ATTACH ref-duckdb (READ_ONLY) and copy:
	##       sequences  (ref `kind=node` rows -- this brings in `seq` so all.fasta is rebuildable)
	##       self_scores
	##       summary_tsv  (already bidirectional from the previous run)
	##   - append input sequences (ord = max_ref_ord + 1..N) via IngestSequences
	##   - re-emit cat/all/all.fasta from ref-duckdb + input_lab2seq (ref first, input after)
	##   - also write cat/all/input.fasta (used by mmseqs branch of 01-2 to build an input subset DB)
	##   - generate per-(input-node, split) tblastx batches when Mode == 'tblastx'
	##   - build the engine DB (makeblastdb / mmseqs createdb) on all.fasta
	if Resume && StepDone.call(t.name)
		PrintStatus.call(args.step, NumStep, "SKIP (resume)", t)
		next
	end
	PrintStatus.call(args.step, NumStep, "START", t)
	jdir     = "#{Odir}/batch/#{t.name.split('.')[0]}"; mkdir_p jdir
	odir     = "#{Odir}/cat/all"; mkdir_p odir
	log      = "#{odir}/all.fasta.makeblastdb.log"
	fa       = "#{odir}/all.fasta"
	input_fa = "#{odir}/input.fasta"
	outs     = []
	puts ""

	## --- 1. validate ref-duckdb compatibility ------------------------------
	ref_meta = {}
	DuckDBUtil.query_each("SELECT key, value FROM run_metadata", db: RefDuckDB) do |line|
		k, v = line.split("\t", 2)
		ref_meta[k] = v
	end
	raise "ref-duckdb has no run_metadata (not a ViPTreeGen run.duckdb?): #{RefDuckDB}" if ref_meta.empty?
	sv = ref_meta["schema_version"]
	if sv.nil? || Gem::Version.new(sv) < Gem::Version.new("2.0")
		raise "ref-duckdb schema_version=#{sv.inspect} is too old; need >= 2.0. " \
		      "Re-generate the reference with ViPTreeGen v2.0+."
	end
	unless ref_meta["mode"] == "normal"
		raise "ref-duckdb was generated in mode=#{ref_meta['mode'].inspect}; only 'normal' refs are supported."
	end
	puts "### ref-duckdb: schema_version=#{sv}, search_mode=#{ref_meta['search_mode']}, " \
	     "matrix=#{ref_meta['matrix']}, evalue=#{ref_meta['evalue']}"
	if ref_meta["search_mode"] != Mode
		puts "### [!] ref-duckdb was built with --mode #{ref_meta['search_mode']} but this run uses --mode #{Mode}. " \
		     "Cross-engine SG values differ slightly; topology should still be preserved."
	end

	## --- 2. parse input fasta ----------------------------------------------
	input_fasta_str = IO.read(Fin)
	input_entries   = ParseFastaEntries.call(input_fasta_str, fasname: "input FASTA file", min_seqs: 1)
	input_labs      = input_entries.map(&:first)
	input_lab2seq   = input_entries.to_h

	## --- 3. read ref ids + max ord from ref-duckdb -------------------------
	ref_ids = {}
	DuckDBUtil.query_each("SELECT seq_id FROM sequences WHERE kind = 'node'", db: RefDuckDB) do |line|
		ref_ids[line.strip] = true
	end
	max_ref_ord = -1
	DuckDBUtil.query_each("SELECT COALESCE(MAX(ord), -1) FROM sequences WHERE kind = 'node'", db: RefDuckDB) do |line|
		max_ref_ord = line.strip.to_i
	end

	## --- 4. detect input × ref ID collisions -------------------------------
	collisions = input_labs.select { |lab| ref_ids[lab] }
	unless collisions.empty?
		raise "input fasta contains seq_id(s) already present in ref-duckdb: " \
		      "#{collisions.first(5).join(', ')}#{collisions.size > 5 ? ' ...' : ''} (#{collisions.size} total).\n" \
		      "Rename or remove the conflicting input sequences."
	end

	## --- 5. initialize new run.duckdb --------------------------------------
	InitDuckDB.call(mode: "ref")

	## --- 6. ATTACH ref-duckdb and copy ref data ----------------------------
	attach_lit = DuckDBUtil.sql_str(RefDuckDB)
	DuckDBUtil.exec_sql(<<~SQL)
		ATTACH #{attach_lit} AS ref_db (READ_ONLY);
		INSERT INTO sequences (seq_id, kind, length, ord, seq)
		  SELECT seq_id, kind, length, ord, seq FROM ref_db.sequences WHERE kind = 'node';
		INSERT INTO self_scores SELECT seq_id, self_score FROM ref_db.self_scores;
		INSERT INTO summary_tsv SELECT * FROM ref_db.summary_tsv;
		DETACH ref_db;
		INSERT INTO run_metadata VALUES ('ref_count', '#{ref_ids.size}');
	SQL

	## --- 7. append input sequences (ord starts after ref) ------------------
	IngestSequences.call(input_labs, input_lab2seq, kind: "node", ord_offset: max_ref_ord + 1)

	## --- 8. emit cat/all/all.fasta, all.len, input.fasta -------------------
	## (ref first, in their original ord; then input)
	open(fa, "w") do |fall|
		open(Flen, "w") do |flen|
			open(input_fa, "w") do |finp|
				DuckDBUtil.query_each(
					"SELECT seq_id, length, seq FROM sequences WHERE kind = 'node' ORDER BY ord",
					db: RefDuckDB,
				) do |line|
					seq_id, length, seq = line.split("\t", 3)
					fall.puts [">"+seq_id, seq.scan(/.{1,70}/)]
					flen.puts [seq_id, length].join("\t")
				end
				input_entries.each do |lab, seq|
					fall.puts [">"+lab, seq.scan(/.{1,70}/)]
					flen.puts [lab, seq.size].join("\t")
					finp.puts [">"+lab, seq.scan(/.{1,70}/)]
				end
			end
		end
	end

	## --- 9. per-(input-node, split) directories + batches -----------------
	## ref nodes get no per-node dir; ref-vs-ref came from ref-duckdb.
	input_entries.each do |lab, seq|
		len   = seq.size
		n0dir = "#{Odir}/node/#{lab}/seq";   mkdir_p n0dir
		n1dir = "#{Odir}/node/#{lab}/blast"; mkdir_p n1dir

		fque = "#{n0dir}/#{lab}.fasta"
		open(fque, "w") { |f| f.puts [">"+lab, seq.scan(/.{1,70}/)] }

		next if UsesMmseqs   # mmseqs runs one global search in 01-2; no per-node batches

		if len > Cutlen
			n0dir_s = "#{n0dir}/split"; mkdir_p n0dir_s
			idx = 0
			while seq && seq.size > 0
				idx += 1
				subseq, seq = seq[0, Cutlen], seq[Cutlen..-1]
				fspt = "#{n0dir_s}/#{idx}.fasta"
				open(fspt, "w") { |f| f.puts [">#{idx}", subseq.scan(/.{1,70}/)] }
				n1dir_s = "#{n1dir}/split/#{idx}"; mkdir_p n1dir_s
				_out = "#{n1dir_s}/tblastx.out"
				_log = "#{n1dir_s}/tblastx.log"
				outs << SearchCmd.call(query: fspt, db: DB, out: _out, log: _log)
			end
		else
			_out = "#{n1dir}/tblastx.out"
			_log = "#{n1dir}/tblastx.log"
			outs << SearchCmd.call(query: fque, db: DB, out: _out, log: _log)
		end
	end

	## --- 10. write batch file (input-only) ---------------------------------
	WriteBatch.call(outs, jdir, t)

	## --- 11. build engine DB on combined all.fasta -------------------------
	sh BuildDBCmd.call(mode: Mode, fasta: Fa, db_path: DB, log: log, create_index: UsesMmseqs)
	DuckDBUtil.ingest_log(log, source: (UsesMmseqs ? "mmseqs:all.createdb" : "makeblastdb"))
	MarkStepDone.call(t.name)
end

desc "01-1.prep_for_tblastx"
task "01-1.prep_for_tblastx", ["step"] do |t, args|
	if Resume && StepDone.call(t.name)
		PrintStatus.call(args.step, NumStep, "SKIP (resume)", t)
		next
	end
	PrintStatus.call(args.step, NumStep, "START", t)
	jdir     = "#{Odir}/batch/#{t.name.split(".")[0]}"; mkdir_p jdir
	odir     = "#{Odir}/cat/all"; mkdir_p odir
	log      = "#{odir}/all.fasta.makeblastdb.log"
	fa       = "#{odir}/all.fasta"
	lab2seq  = {}
	outs     = []
	puts ""

	## initialize DuckDB and record run metadata
	InitDuckDB.call(mode: "normal")

	## validate and make copy of input fasta
	open(Flen, "w"){ |flen|
		open(fa, "w"){ |fall|
			fasta_str = IO.read(Fin)
			entries   = ParseFastaEntries.call(fasta_str, fasname: "input FASTA file", min_seqs: 3)

			entries.each{ |lab, seq|
				len = seq.size

				# store sequence name
				lab2seq[lab] = seq

				# make cat output
				flen.puts [lab, len]*"\t"
				fall.puts [">"+lab, seq.scan(/.{1,70}/)]

				# make node output
				# split fasta and make tblastx/blastn job
				n0dir = "#{Odir}/node/#{lab}/seq";   mkdir_p n0dir
				n1dir = "#{Odir}/node/#{lab}/blast"; mkdir_p n1dir

				fque = "#{n0dir}/#{lab}.fasta"
				open(fque, "w"){ |_fque|
					_fque.puts [">"+lab, seq.scan(/.{1,70}/)]
				}

				## mmseqs runs one global search in 01-2 -- skip per-(node, split) batch generation.
				## m8_split.rb in 01-2 will create blast/tblastx.out files on demand.
				next if UsesMmseqs

				if len > Cutlen
					n0dir_s = "#{n0dir}/split"; mkdir_p n0dir_s

					idx = 0
					while seq and seq.size > 0
						idx += 1
						subseq, seq = seq[0, Cutlen], seq[Cutlen..-1]
						fspt = "#{n0dir_s}/#{idx}.fasta"
						open(fspt, "w"){ |_fspt|
							_fspt.puts [">#{idx}", subseq.scan(/.{1,70}/)]
						}

						n1dir_s = "#{n1dir}/split/#{idx}"; mkdir_p n1dir_s
						_out  = "#{n1dir_s}/tblastx.out"
						_log  = "#{n1dir_s}/tblastx.log"
						outs << SearchCmd.call(query: fspt, db: DB, out: _out, log: _log)
					end
				else
					_out = "#{n1dir}/tblastx.out"
					_log = "#{n1dir}/tblastx.log"
					outs << SearchCmd.call(query: fque, db: DB, out: _out, log: _log)
				end
			}
		}
	}

	## write batch file
	WriteBatch.call(outs, jdir, t)

	## build target DB for all.fasta (createindex for mmseqs to amortize per-query lookups)
	sh BuildDBCmd.call(mode: Mode, fasta: Fa, db_path: DB, log: log, create_index: (UsesMmseqs))

	## populate DuckDB: sequences + makeblastdb.log (file kept for backward compat in P0)
	IngestSequences.call(lab2seq.keys, lab2seq, kind: "node")
	DuckDBUtil.ingest_log(log, source: (UsesMmseqs ? "mmseqs:all.createdb" : "makeblastdb"))
	MarkStepDone.call(t.name)
end
desc "01-2.tblastx"
task "01-2.tblastx", ["step"] do |t, args|
	## Run the all-vs-all search using the selected engine (see `Mode`).
	## tblastx mode:   RunBatch over the per-(node,split) batch files written by 01-1.
	## mmseqs mode:    1 global `mmseqs search` (using its internal multi-threading),
	##                 then `convertalis` to BLAST-tab, then `m8_split.rb` to lay out
	##                 the result as `<kind>/<qid>/blast/tblastx.out` files so that
	##                 the Rust binary in step 02-3 reads them unchanged.
	if Resume && StepDone.call(t.name)
		PrintStatus.call(args.step, NumStep, "SKIP (resume)", t)
		next
	end
	PrintStatus.call(args.step, NumStep, "START", t)

	case Mode
	when "tblastx", "blastn"
		jdir = "#{Odir}/batch/01-1" # output from 01-1
		RunBatch.call(jdir, Qname, Nthreads, Mem, Wtime, Ncpus)

	when "mmseqs-tblastx"
		## mmseqs --search-type 4 = translated 6-frame protein alignment. Defaults BLOSUM62,
		## but ViPTreeGen targets divergent viral genomes for which BLOSUM45 is more sensitive,
		## so we ship blosum45.out (mmseqs2 format, freely redistributable) under data/ and pass
		## it via --sub-mat. This also keeps tblastx and mmseqs-tblastx on the same substitution
		## matrix so their SG scores agree tightly. (mmseqs-blastn was evaluated and removed --
		## see doc/mmseqs-blastn.md.)
		require "fileutils"
		threads = Ncpus != "" ? Ncpus : "1"
		tmpdir  = "#{Odir}/mmseqs_tmp"; FileUtils.mkdir_p(tmpdir)
		fmt = '"query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"'
		m8_split = "#{File.dirname(__FILE__)}/script/m8_split.rb"
		sub_mat_path = "#{File.dirname(__FILE__)}/data/#{Matrix.downcase}.out"
		unless File.file?(sub_mat_path)
			raise "\e[1;31mError:\e[0m mmseqs substitution matrix not found: #{sub_mat_path}. " \
			      "Bundled matrices live under data/ in the repo; ensure --matrix matches one of them."
		end
		## `--split-memory-limit X` caps the peak RAM mmseqs uses for the target index; if the
		## index would exceed X, mmseqs auto-splits the search into multiple passes. Default 12G.
		search_opts = "--search-type 4 -e #{Evalue} --max-seqs #{Max_target_seqs} " \
		              "--threads #{threads} --split-memory-limit #{MmseqsSplitMem} " \
		              "--sub-mat #{sub_mat_path} -v 1"

		## per-step mmseqs logs land in cat/all/ so FinalLogCleanup can ingest them into the `logs` table.
		ldir = "#{Odir}/cat/all"

		if RefMode
			## with-reference mode: build a query DB from input.fasta only, then run
			## `mmseqs search input_db all_db` so the result HSPs cover input×input and
			## input×ref pairs (ref×ref came from ref-duckdb already, via 01-1.ref.prep).
			input_fa = "#{ldir}/input.fasta"
			input_db = "#{Fa}.input.mmdb"
			res = "#{tmpdir}/result"
			m8  = "#{tmpdir}/result.m8"
			sh BuildDBCmd.call(mode: Mode, fasta: input_fa, db_path: input_db,
			                   log: "#{input_fa}.createdb.log", create_index: false)
			sh "mmseqs search #{input_db} #{DB} #{res} #{tmpdir}/s #{search_opts} >#{ldir}/mmseqs_search.log 2>&1"
			sh "mmseqs convertalis #{input_db} #{DB} #{res} #{m8} --format-output #{fmt} >#{ldir}/mmseqs_convertalis.log 2>&1"
			sh "ruby #{m8_split} #{m8} #{Odir}/node"
		elsif File.exist?(Flen_i)
			## 2D mode: query (node) vs all, and input vs input self-blast.
			query_db = "#{Fa}.query.mmdb"
			input_db = "#{Fa}.input.mmdb"
			res_node  = "#{tmpdir}/result.node"
			res_input = "#{tmpdir}/result.input"
			m8_node   = "#{tmpdir}/result.node.m8"
			m8_input  = "#{tmpdir}/result.input.m8"

			sh "mmseqs search #{query_db} #{DB} #{res_node} #{tmpdir}/n #{search_opts} >#{ldir}/mmseqs_search.node.log 2>&1"
			sh "mmseqs convertalis #{query_db} #{DB} #{res_node} #{m8_node} --format-output #{fmt} >#{ldir}/mmseqs_convertalis.node.log 2>&1"
			sh "mmseqs search #{input_db} #{input_db} #{res_input} #{tmpdir}/i #{search_opts} >#{ldir}/mmseqs_search.input.log 2>&1"
			sh "mmseqs convertalis #{input_db} #{input_db} #{res_input} #{m8_input} --format-output #{fmt} >#{ldir}/mmseqs_convertalis.input.log 2>&1"
			sh "ruby #{m8_split} #{m8_node}  #{Odir}/node"
			sh "ruby #{m8_split} #{m8_input} #{Odir}/input"
		else
			## normal mode: all-vs-all in one search.
			res = "#{tmpdir}/result"
			m8  = "#{tmpdir}/result.m8"
			sh "mmseqs search #{DB} #{DB} #{res} #{tmpdir}/s #{search_opts} >#{ldir}/mmseqs_search.log 2>&1"
			sh "mmseqs convertalis #{DB} #{DB} #{res} #{m8} --format-output #{fmt} >#{ldir}/mmseqs_convertalis.log 2>&1"
			sh "ruby #{m8_split} #{m8} #{Odir}/node"
		end
		sh "rm -rf #{tmpdir}"
	end
	MarkStepDone.call(t.name)
end
# }}} tasks 01


# {{{ tasks 02
desc "02-3.make_summary_pre"
task "02-3.make_summary_pre", ["step"] do |t, args|
	## Reads per-node tblastx.out files (node/<qid>/blast/{split/<i>/,}tblastx.out)
	## and populates the summary_pre table via the Rust binary
	## viptreegen-summary-pre, which fuses the legacy 01-3 (cat+rename+shift),
	## 02-1 (idt/alen filter), 02-2 (BED conversion), and 02-3 (interval merge)
	## into one parallel single-process step.
	if Resume && StepDone.call(t.name)
		PrintStatus.call(args.step, NumStep, "SKIP (resume)", t)
		next
	end
	PrintStatus.call(args.step, NumStep, "START", t)
	tdir    = "#{Odir}/tmp"; mkdir_p tdir
	bin     = RustBin.call
	threads = (Ncpus != "" ? Ncpus : "1")
	idt     = (ENV["idt"]   || "30")
	alen    = (ENV["aalen"] || "30")

	kinds = File.exist?(Flen_i) ? %w|node input| : %w|node|
	sql = +"BEGIN;\nDELETE FROM summary_pre;\n"
	kinds.each do |kind|
		out_tsv = "#{tdir}/02-3.#{kind}.tsv"
		sh "#{bin} --outdir #{Odir} --kind #{kind} --cutlen #{Cutlen} " \
		   "--idt #{idt} --alen #{alen} --threads #{threads} --output #{out_tsv}"
		next if File.zero?(out_tsv)
		sql << "COPY summary_pre FROM '#{out_tsv}' (FORMAT CSV, DELIMITER E'\\t', HEADER false);\n"
	end
	sql << "CREATE INDEX IF NOT EXISTS idx_sp_kind_node ON summary_pre(kind, node);\n"
	sql << "COMMIT;\n"
	DuckDBUtil.exec_sql(sql)

	## cleanup: tblastx.out files (now consumed into summary_pre) + tmp TSVs
	sh "find #{Odir}/node #{Odir}/input -type f -name 'tblastx.out' -delete 2>/dev/null || true; rm -f #{tdir}/02-3.*.tsv"
	MarkStepDone.call(t.name)
end
desc "02-4.make_self_tblastx"
task "02-4.make_self_tblastx", ["step"] do |t, args|
	if Resume && StepDone.call(t.name)
		PrintStatus.call(args.step, NumStep, "SKIP (resume)", t)
		next
	end
	PrintStatus.call(args.step, NumStep, "START", t)

	## populate self_scores via pure SQL: self-hit is where que == sub.
	## In 2D mode the same seq_id appears under both kind='node' and kind='input';
	## we prefer the 'input' value (matches legacy iteration order: node first, input last wins).
	##
	## In ref mode self_scores already contains ref entries imported from --ref-duckdb
	## (01-1.ref.prep). We only want to ADD input self_scores. summary_pre holds only
	## input × * HSPs in ref mode, so the INSERT naturally scopes to input ids.
	delete_clause = RefMode ? "" : "DELETE FROM self_scores;\n"
	DuckDBUtil.exec_sql <<~SQL
		BEGIN;
		#{delete_clause}INSERT INTO self_scores (seq_id, self_score)
			SELECT que, (que_score + sub_score) * 0.5 FROM (
				SELECT que, que_score, sub_score,
				       ROW_NUMBER() OVER (PARTITION BY que ORDER BY CASE kind WHEN 'input' THEN 0 ELSE 1 END) AS rn
				FROM summary_pre WHERE que = sub
			) WHERE rn = 1;
		COMMIT;
	SQL
	MarkStepDone.call(t.name)
end
desc "02-5.2D.make_summary_tsv" ### 2D mode
task "02-5.2D.make_summary_tsv", ["step"] do |t, args|
	if Resume && StepDone.call(t.name)
		PrintStatus.call(args.step, NumStep, "SKIP (resume)", t)
		next
	end
	PrintStatus.call(args.step, NumStep, "START", t)
	script   = "#{File.dirname(__FILE__)}/script/#{t.name}.rb"
	tdir     = "#{Odir}/tmp"; mkdir_p tdir
	out_tsv  = "#{tdir}/02-5.all.tsv"
	sh "ruby #{script} #{out_tsv}"

	DuckDBUtil.exec_sql <<~SQL
		BEGIN;
		DELETE FROM summary_tsv;
		COPY summary_tsv FROM '#{out_tsv}' (FORMAT CSV, DELIMITER E'\\t', HEADER false);
		CREATE INDEX IF NOT EXISTS idx_stsv_node ON summary_tsv(node);
		COMMIT;
	SQL
	sh "rm -f #{out_tsv}"
	MarkStepDone.call(t.name)
end
desc "02-5.make_summary_tsv" ### normal mode (also used by ref mode)
task "02-5.make_summary_tsv", ["step"] do |t, args|
	if Resume && StepDone.call(t.name)
		PrintStatus.call(args.step, NumStep, "SKIP (resume)", t)
		next
	end
	PrintStatus.call(args.step, NumStep, "START", t)
	script   = "#{File.dirname(__FILE__)}/script/#{t.name}.rb"
	tdir     = "#{Odir}/tmp"; mkdir_p tdir
	out_tsv  = "#{tdir}/02-5.all.tsv"
	## In ref mode the script needs to know which seq_ids are refs so that input×ref HSPs
	## (only one direction exists in summary_pre) are not discarded by the legacy
	## shorter→longer dedup filter. Pass ref count via env so the script can build the set.
	ENV["ref_mode"] = RefMode ? "true" : ""
	sh "ruby #{script} #{out_tsv}"

	## In ref mode summary_tsv already contains ref×ref rows (imported from --ref-duckdb).
	## summary_pre only holds input × * HSPs, so the new rows append cleanly without
	## colliding with the imported ref×ref rows. Skip DELETE to preserve them.
	delete_clause = RefMode ? "" : "DELETE FROM summary_tsv;\n"
	DuckDBUtil.exec_sql <<~SQL
		BEGIN;
		#{delete_clause}COPY summary_tsv FROM '#{out_tsv}' (FORMAT CSV, DELIMITER E'\\t', HEADER false);
		CREATE INDEX IF NOT EXISTS idx_stsv_node ON summary_tsv(node);
		COMMIT;
	SQL
	sh "rm -f #{out_tsv}"
	MarkStepDone.call(t.name)
end
# }}} tasks 02


# {{{ tasks 03
desc "03-1.2D.make_matrix"
task "03-1.2D.make_matrix", ["step"] do |t, args|
	if Resume && StepDone.call(t.name)
		PrintStatus.call(args.step, NumStep, "SKIP (resume)", t)
		next
	end
	PrintStatus.call(args.step, NumStep, "START", t)
	rdir = "#{Odir}/result"; mkdir_p rdir
	qids = []; DuckDBUtil.query_each("SELECT seq_id FROM sequences WHERE kind = 'node' ORDER BY ord") { |l| qids << l.strip }
	tids = []; DuckDBUtil.query_each("SELECT seq_id FROM sequences WHERE kind = 'input' ORDER BY ord") { |l| tids << l.strip }

	similarities = Hash.new{ |h, i| h[i] = Hash.new(0.0) }
	DuckDBUtil.query_each("SELECT node, target, sg FROM summary_tsv"){ |l|
		id1, id2, sim = l.split("\t")
		similarities[id1][id2] = sim.to_f
		similarities[id2][id1] = sim.to_f
	}

	outs = []
	qids.each{ |id1|
		out = [id1]
		tids.each{ |id2|
			s = id1 == id2 ? 1 : similarities[id1][id2]
			s = s == 0 ? "0" : (s == 1 ? "1" : "%.4f" % s)
			out << s
		}
		outs << out
	}
	open("#{rdir}/2D.sim.matrix", "w"){ |fout|
		fout.puts ["", tids]*"\t"
		outs.each{ |out|
			fout.puts out*"\t"
		}
	}
	open("#{rdir}/top10.sim.list", "w"){ |fout|
		outs.each{ |out|
			query, *scores = out
			score2tids = Hash.new{ |h, i| h[i] = [] }
			tids.zip(scores){ |tid, score|
				score2tids[score] << tid
			}
			a = []
			score2tids.sort_by{ |score, tids| -score.to_f }.each{ |score, tids|
				tids.each{ |tid|
					break if a.size >= 10
					a << "#{tid}:#{score}"
				}
			}
			fout.puts [query, a]*"\t"
		}
	}

	FinalLogCleanup.call
	MarkStepDone.call(t.name)
end
desc "03-1.make_matrix"
task "03-1.make_matrix", ["step"] do |t, args|
	if Resume && StepDone.call(t.name)
		PrintStatus.call(args.step, NumStep, "SKIP (resume)", t)
		next
	end
	PrintStatus.call(args.step, NumStep, "START", t)
	rdir = "#{Odir}/result"; mkdir_p rdir
	ids  = []; DuckDBUtil.query_each("SELECT seq_id FROM sequences WHERE kind = 'node' ORDER BY ord") { |l| ids << l.strip }

	similarities = Hash.new{ |h, i| h[i] = Hash.new(0.0) }
	fout1   = open("#{rdir}/all.sim.matrix", "w")
	fout2   = open("#{rdir}/all.dist.matrix", "w")

	DuckDBUtil.query_each("SELECT node, target, sg FROM summary_tsv"){ |l|
		id1, id2, sim = l.split("\t")
		similarities[id1][id2] = sim.to_f
	}

	[fout1, fout2].each{ |fout| fout.puts ["", ids]*"\t" }
	ids.each{ |id1|
		out1 = [id1]
		out2 = [id1]
		ids.each{ |id2|
			s = similarities[id1][id2]
			s = s == 0 ? "0" : (s == 1 ? "1" : "%.4f" % s)
			d = 1 - similarities[id1][id2]
			d = d == 0 ? "0" : (d == 1 ? "1" : "%.4f" % d)
			out1 << s
			out2 << d
		}
		fout1.puts out1*"\t"
		fout2.puts out2*"\t"
	}
	[fout1, fout2].each{ |fout| fout.close }

	FinalLogCleanup.call
	MarkStepDone.call(t.name)
end
desc "03-2.matrix_to_nj"
task "03-2.matrix_to_nj", ["step"] do |t, args|
	if Resume && StepDone.call(t.name)
		PrintStatus.call(args.step, NumStep, "SKIP (resume)", t)
		next
	end
	PrintStatus.call(args.step, NumStep, "START", t)
	rdir     = "#{Odir}/result"
	flag     = "dist"
	script   = "#{File.dirname(__FILE__)}/script/#{t.name}.R"
	fin      = "#{rdir}/all.#{flag}.matrix"
	foutpref = "#{rdir}/all.#{TreeMethod}"
	flog     = "#{rdir}/all.#{TreeMethod}.makelog"
	sh "LANG=C Rscript --quiet --no-save --no-restore #{script} #{fin} #{foutpref} #{flag} #{TreeMethod} >#{flog} 2>&1"
	MarkStepDone.call(t.name)
end
# }}} tasks 03
