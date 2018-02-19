#
#  02-5.make_summary_tsv.rb / ViPTreeGen - a command line tool for viral proteomic tree generation
#
#    Copyright: 2017 (C) Yosuke Nishimura (yosuke@kuicr.kyoto-u.ac.jp)
#    License: MIT license
#    Initial version: 2017-02-07
#

flen, dir, type = ARGV
header = %w|ID length score SG %.mean.idt %.len|

id2info  = {}
id2order = {}
IO.readlines(flen).each{ |l| 
	id, len = l.chomp.split("\t")
	id2info[id]  = [id, len]
	id2order[id] = [len.to_i, id]
}

fins   = Dir["#{dir}/node/*/blast/#{type}.summary.pre"]
id2out = Hash.new{ |h, i| h[i] = [] }
size   = fins.size
fins.each.with_index(1){ |fin, idx|
	#              1               2    3    4          5          6        7        8               9              10
	# que_score_rank  sub_score_rank  que  sub  que_score  sub_score  que_len  sub_len  que_len_in_hit  sub_len_in_hit 
	# 						  11                12               13               14
	# %mean_idt_of_que  %mean_idt_of_sub  %que_len_in_hit  %sub_len_in_hit
	IO.readlines(fin)[1..-1].each{ |l|
		a = l.chomp.split("\t")[2..-1]
		que, sub = a[0..1]

		# USE ONLY HSP of [SHORTER => LONGER]
		# next if id2order[que] > id2order[sub]
		next if id2order[que] <=> id2order[sub] == 1 ### compare [length (as integer), id (in dictionary sort)]

		id2out[que] << a
		id2out[sub] << a if que != sub
	}
	print "#{idx} / #{size} finished" if idx % 1000 == 0
}


# self score for SG
id2self = {}
fin1    = "#{dir}/cat/blast/#{type}.self.sum"
IO.readlines(fin1).each{ |l|
	id, score   = l.chomp.split("\t")
	id2self[id] = score.to_f
}

fins.each{ |fin|
	id = fin.split("/")[2]
	outs = []
	id2out[id].each{ |a|
		a = a.values_at(0..3, 8..11) # q, s, q_scr, s_scr, q_%idt, s_%idt, q_%len, s_%len
		query_is_id = a[0] == id
		targ  = query_is_id ? a[1] : a[0]
		p_idt = query_is_id ? a[4] : a[5]
		p_len = query_is_id ? a[6] : a[7]
		score = (a[2].to_i + a[3].to_i) * 0.5
		s_scr = query_is_id ? id2self[id] : id2self[targ] # self score of shorter sequence
		sg = "%.4f" % [score / s_scr, 1].min              # max of SG should be 1
		outs << [targ, score, sg, p_idt, p_len]
	}

	open("#{dir}/node/#{id}/blast/#{type}.summary.tsv", "w"){ |fout|
		fout.puts header*"\t"
		outs.sort_by{ |a| -a[2].to_f }.each{ |a| fout.puts [id2info[a[0]], a[1..-1]]*"\t" }
	}
}
