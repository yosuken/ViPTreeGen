#
#  01-3.cat_and_rename_split_tblastx.rb / ViPTreeGen - a command line tool for viral proteomic tree generation
#
#    Copyright: 2017 (C) Yosuke Nishimura (yosuke@kuicr.kyoto-u.ac.jp)
#    License: MIT license
#    Initial version: 2017-02-07
#

node, program, fa, n1dir, cutlen = ARGV
cutlen = cutlen.to_i

out   = []
Dir["#{n1dir}/split/*/#{program}.out"].sort_by{ |fin| fin.split("/")[-2].to_i }.each{ |fin|
	idx = fin.split("/")[-2].to_i - 1
	open("#{fin}2", "w"){ |fout|
		IO.readlines(fin).each{ |l| 
			a = l.chomp.split("\t")
			a[0] = node # change query name
			(6..7).each{ |i| a[i] = a[i].to_i + cutlen * idx } # q.start, q.end is moved by +(cutlen * idx)nt.
			out << a
			fout.puts a*"\t"
		}
	}
}

# order = IO.readlines(flen).each_with_index.inject({}){ |h, (l, idx)| h[l.split("\t")[0]] = idx; h }
order = {} ###
IO.read(fa).split(/^>/)[1..-1].each_with_index{ |ent, idx|
	lab, *seq = ent.split("\n") 
	order[lab] = idx
}

open("#{n1dir}/#{program}.out", "w"){ |fout|
	out.sort_by{ |a| [order[a[1]], a[10].to_f] }.each{ |a| fout.puts a*"\t" }
}
