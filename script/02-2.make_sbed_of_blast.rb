#!/usr/bin/env ruby
#
#  02-2.make_sbed_of_blast.rb / ViPTreeGen - a command line tool for viral proteomic tree generation
#
#    Copyright: 2017 (C) Yosuke Nishimura (yosuke@kuicr.kyoto-u.ac.jp)
#    License: MIT license
#    Initial version: 2017-02-07
#

# Convert blast tabular to BED6 file
# blast tabular: qID, sID, %identity, alignLen, mismatches, gapOpens, q.start, q.end, s.start, s.end, evalue, bitScore
# BED6         : sID, s.start, s.end, qID:q.start-q.end, identity, strand
# [identity] 99.1% => 991
# [positions] 0-based start positions; s.start < s.end; q.start < q.end

# e.g.
# input: Q1 S1 93.41 744 44 5 938 198 3231 2490 0.0 1098
# S1 2489 3231 Q1:197-938 934 +

fin, fout = ARGV
fout = open(fout, "w")
out = []

IO.readlines(fin).each.with_index(1){ |l, idx|
	next if l =~ /^#/
	a = l.chomp.split(/\t/)
	raise "invalid input file: # of column is detected as #{a.size} in line #{idx}; #{l.inspect}" if a.size != 12

	query,     subject   = a[0],      a[1]
	query_e,   query_s   = a[7].to_i, a[6].to_i
	subject_e, subject_s = a[9].to_i, a[8].to_i
	strand               = (subject_e - subject_s) * (query_e - query_s) > 0 ? "+" : "-"

	subject_e, subject_s = subject_s,   subject_e if subject_s > subject_e
	query_e,   query_s   = query_s,     query_e   if query_s > query_e

	# blast(1-based): [1, 100] => bed(0-based): [0, 100)
	query_s,   subject_s = query_s - 1, subject_s - 1 

	# percent_identity * 10 (int, 0 ~ 1000)
	score = (a[2].to_f * 10).round 

	query_label = [query, ":", query_s, "-", query_e]*""
	o = [subject, subject_s, subject_e, query_label, score, strand]

	# additional 7th column: bitScore
	bit_score = a[11]
	o += [bit_score]

	out << o
}
out.sort_by{ |o| [o[0], o[1].to_i, o[2].to_i, o[3]] }.each{ |o| fout.puts o.join("\t") }
fout.close
