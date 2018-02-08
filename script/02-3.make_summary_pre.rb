#
#  02-3.make_summary_pre.rb / ViPTreeGen - a command line tool for viral proteomic tree generation
#
#    Copyright: 2017 (C) Yosuke Nishimura (yosuke@kuicr.kyoto-u.ac.jp)
#    License: MIT license
#    Initial version: 2017-02-07
#

# fin, fout, node, flen = ARGV
n1dir, fout, node, flen, filename = ARGV
fout = open(fout, "w")

id2len = {}
IO.readlines(flen).each{ |l| id, len = l.chomp.split("\t"); id2len[id] = len }

lab2sub_pos_idt = Hash.new{ |h, i| h[i] = {} } # lab2sub_pos_idt[lab] = { sub_range => identity }
lab2sub_pos_scr = Hash.new{ |h, i| h[i] = {} } # lab2sub_pos_scr[lab] = { sub_range => score    }
lab2que_pos_idt = Hash.new{ |h, i| h[i] = {} } # lab2que_pos_idt[lab] = { que_range => identity }
lab2que_pos_scr = Hash.new{ |h, i| h[i] = {} } # lab2que_pos_scr[lab] = { que_range => score    }

fins = File.exist?("#{n1dir}/#{filename}") ? ["#{n1dir}/#{filename}"] : Dir["#{n1dir}/split/*/#{filename}"]
fins.each{ |fin|
	IO.readlines(fin).each{ |l|
		lab, start, stop, que, identity, strand, score = l.chomp.split("\t")

		# subject start, stop
		start, stop = start.to_i + 1, stop.to_i # 0-based => 1-based
		raise if start > stop
		lab2sub_pos_idt[lab][start..stop] = identity.to_f
		lab2sub_pos_scr[lab][start..stop] = score.to_f / (stop - start + 1)

		# query  start, stop
		start, stop = que.split(":")[1].split("-")
		start, stop = start.to_i + 1, stop.to_i # 0-based => 1-based
		raise if start > stop
		lab2que_pos_idt[lab][start..stop] = identity.to_f
		lab2que_pos_scr[lab][start..stop] = score.to_f / (stop - start + 1)
	}
}
raise if lab2sub_pos_idt.keys != lab2que_pos_idt.keys
raise if lab2sub_pos_idt.keys != lab2sub_pos_scr.keys
raise if lab2que_pos_scr.keys != lab2sub_pos_scr.keys

## strategy
# [in range]
# 1.......9  
#  2..5       
#    4...8  
#    4......B
# [points]
# 12 45  89 B

lab2sub_point = Hash.new{ |h, i| h[i] = Hash.new{ |g, j| g[j] = [] } } # lab2sub_point[lab] = { point => range }
lab2que_point = Hash.new{ |h, i| h[i] = Hash.new{ |g, j| g[j] = [] } } # lab2que_point[lab] = { point => range }

# subject
lab2sub_pos_idt.each{ |lab, pos_h|
	ranges = pos_h.keys
	points = ranges.map{ |r| [r.begin, r.end] }.flatten.sort
	points.each{ |point|
		ranges.each{ |range|
			lab2sub_point[lab][point] << range if range.begin <= point and range.end >= point
		}
	}
}
# query
lab2que_pos_idt.each{ |lab, pos_h|
	ranges = pos_h.keys
	points = ranges.map{ |r| [r.begin, r.end] }.flatten.sort
	points.each{ |point|
		ranges.each{ |range|
			lab2que_point[lab][point] << range if range.begin <= point and range.end >= point
		}
	}
}

lab2sub_pos_idt2 = Hash.new{ |h, i| h[i] = {} } # lab2sub_pos_idt2[lab] = { range => max_identity } (splited by points)
lab2sub_pos_scr2 = Hash.new{ |h, i| h[i] = {} } # lab2sub_pos_scr2[lab] = { range => max_score    } (splited by points)
lab2que_pos_idt2 = Hash.new{ |h, i| h[i] = {} } # lab2que_pos_idt2[lab] = { range => max_identity } (splited by points)
lab2que_pos_scr2 = Hash.new{ |h, i| h[i] = {} } # lab2que_pos_scr2[lab] = { range => max_score    } (splited by points)

pos_a  = [[lab2sub_pos_idt,  lab2sub_pos_scr],  [lab2que_pos_idt,  lab2que_pos_scr]]
pos2_a = [[lab2sub_pos_idt2, lab2sub_pos_scr2], [lab2que_pos_idt2, lab2que_pos_scr2]]
# subject
lab2sub_point.each{ |lab, sub_point_h|
	[lab2sub_pos_idt, lab2sub_pos_scr].zip([lab2sub_pos_idt2, lab2sub_pos_scr2]){ |lab2sub_pos, lab2sub_pos2|
		points = sub_point_h.keys.sort
		(0..points.size - 1).each{ |i|
			# best for point i 
			i_ranges  = lab2sub_point[lab][points[i]]
			i_max = i_ranges.map{ |range| lab2sub_pos[lab][range] }.max
			lab2sub_pos2[lab][points[i]..points[i]] = i_max

			# best for range of (i, i+1)
			next if i == points.size - 1 # end point
			j = i + 1
			j_ranges  = lab2sub_point[lab][points[j]]
			b_ranges  = i_ranges & j_ranges    # in-between ranges = shared ranges
			next if b_ranges.size == 0         # skip gap
			b_max = b_ranges.map{ |range| lab2sub_pos[lab][range] }.max
			start, stop = points[i] + 1, points[j] - 1

			lab2sub_pos2[lab][start..stop] = b_max if start <= stop
		}
	}
}
# query
lab2que_point.each{ |lab, que_point_h|
	[lab2que_pos_idt, lab2que_pos_scr].zip([lab2que_pos_idt2, lab2que_pos_scr2]){ |lab2que_pos, lab2que_pos2|
		points = que_point_h.keys.sort
		(0..points.size - 1).each{ |i|
			# best for point i 
			i_ranges  = lab2que_point[lab][points[i]]
			i_max = i_ranges.map{ |range| lab2que_pos[lab][range] }.max
			lab2que_pos2[lab][points[i]..points[i]] = i_max

			# best for range of (i, i+1)
			next if i == points.size - 1 # end point
			j = i + 1
			j_ranges  = lab2que_point[lab][points[j]]
			b_ranges  = i_ranges & j_ranges    # in-between ranges = shared ranges
			next if b_ranges.size == 0         # skip gap
			b_max = b_ranges.map{ |range| lab2que_pos[lab][range] }.max
			start, stop = points[i] + 1, points[j] - 1

			lab2que_pos2[lab][start..stop] = b_max if start <= stop
		}
	}
}

out = []
labels = lab2sub_pos_idt.keys
#labels = lab2sub_point.keys
#(p (labels - lab2que_point.keys); p (lab2que_point.keys - labels); raise) if labels != lab2que_point.keys

labels.each{ |lab|
	total_sub_len = 0
	total_sub_idt = 0
	total_sub_scr = 0
	total_que_len = 0
	total_que_idt = 0
	total_que_scr = 0

	sub_idt_h = lab2sub_pos_idt2[lab]
	sub_scr_h = lab2sub_pos_scr2[lab]
	que_idt_h = lab2que_pos_idt2[lab]
	que_scr_h = lab2que_pos_scr2[lab]
	#p sub_idt_h

	sub_set = [sub_idt_h, sub_scr_h, total_sub_len, total_sub_idt, total_sub_scr]
	que_set = [que_idt_h, que_scr_h, total_que_len, total_que_idt, total_que_scr]
	sub_idt_h.each_key{ |range|
		len = range.size
		total_sub_len += len
		total_sub_idt += len * sub_idt_h[range]
		total_sub_scr += len * sub_scr_h[range]
	}
	que_idt_h.each_key{ |range|
		len = range.size
		total_que_len += len
		total_que_idt += len * que_idt_h[range]
		total_que_scr += len * que_scr_h[range]
	}
	#p total_sub_scr
	#p total_sub_len

	que             = node
	sub             = lab
	que_len         = id2len[que]
	sub_len         = id2len[sub]
	mean_que_idt    = "%.1f" % (total_que_idt.to_f / total_que_len / 10) # permil => percent
	mean_sub_idt    = "%.1f" % (total_sub_idt.to_f / total_sub_len / 10) # permil => percent
	total_que_scr   = total_que_scr.round
	total_sub_scr   = total_sub_scr.round
	percent_que_len = "%.1f" % (total_que_len.to_f / que_len.to_f * 100)
	percent_sub_len = "%.1f" % (total_sub_len.to_f / sub_len.to_f * 100)
	out << [que, sub, total_que_scr, total_sub_scr, que_len, sub_len, total_que_len, total_sub_len, mean_que_idt, mean_sub_idt, percent_que_len, percent_sub_len]
}

# header
header  = %w|que_score_rank sub_score_rank que sub|
header += %w|que_score sub_score que_len sub_len que_len_in_hit sub_len_in_hit %mean_idt_of_que %mean_idt_of_sub %que_len_in_hit %sub_len_in_hit|
fout.puts header*"\t"

# sort then write
out = out.sort_by{ |a| -a[3].to_f }.map.with_index(1){ |a, idx| [idx] + a } # sub_score_rank
out = out.sort_by{ |a| -a[3].to_f }.map.with_index(1){ |a, idx| [idx] + a } # que_score_rank
out.each{ |a| fout.puts a*"\t" }
