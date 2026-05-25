#
#  m8_split.rb / ViPTreeGen
#
#    Copyright: 2026 (C) Yosuke Nishimura
#    License: MIT
#
#  Split a single mmseqs convertalis `.m8` (BLAST-tab) file by the first column
#  (query id) into per-qid files at `<outroot>/<qid>/blast/tblastx.out`.
#  mmseqs2 emits hits in query-grouped order (each query's hits come together),
#  so we stream-split with a single open file handle at a time.
#
#  Usage:
#    ruby m8_split.rb <m8_file> <out_root_dir>
#
#  The Rust binary viptreegen-summary-pre then discovers
#  `<outroot>/<qid>/blast/tblastx.out` (no split subdir) and treats them as
#  non-split per-node files.

require "fileutils"

m8, outroot = ARGV
raise "Usage: ruby m8_split.rb <m8_file> <out_root_dir>" unless m8 && outroot

current_qid = nil
current_fh  = nil
nrows = 0
nfiles = 0

IO.foreach(m8) do |line|
  qid = line.split("\t", 2).first
  if qid != current_qid
    current_fh.close if current_fh
    dir = "#{outroot}/#{qid}/blast"
    FileUtils.mkdir_p(dir)
    current_fh = open("#{dir}/tblastx.out", "w")
    current_qid = qid
    nfiles += 1
  end
  current_fh.write(line)
  nrows += 1
end
current_fh.close if current_fh

warn "[m8_split] wrote #{nrows} rows across #{nfiles} qid files under #{outroot}"
