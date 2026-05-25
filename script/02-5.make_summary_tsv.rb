#
#  02-5.make_summary_tsv.rb / ViPTreeGen - a command line tool for viral proteomic tree generation
#
#    Copyright: 2017 (C) Yosuke Nishimura (yosuke@kuicr.kyoto-u.ac.jp)
#    License: MIT license
#    Initial version: 2017-02-07
#
#  Reads summary_pre + self_scores + sequences from DuckDB, computes the per-node
#  (target, length, score, sg, mean_idt, pct_len) summary, and writes a single combined
#  TSV which the rake task COPYs into the `summary_tsv` table.
#
#  Args: out_tsv
#

require_relative "duckdb_util"
require "set"

out_tsv = ARGV[0]

## id2info from all sequences (covers both kinds in case targets span them); id2order same.
id2info  = {}
id2order = {}
DuckDBUtil.query_each("SELECT DISTINCT seq_id, length FROM sequences ORDER BY ord"){ |l|
	id, len = l.split("\t")
	id2info[id]  = [id, len]
	id2order[id] = [len.to_i, id]
}

## sequences iterated for output (only kind='node' owns a summary tsv row group)
node_ids = []
DuckDBUtil.query_each("SELECT seq_id FROM sequences WHERE kind = 'node' ORDER BY ord"){ |l|
	node_ids << l.strip
}

## In ref mode the first `ref_count` node ids (in ord) are reference sequences. Build the set
## so the legacy shorter→longer filter can be skipped for input × ref pairs (only one direction
## of those HSPs exists in summary_pre, so the filter would otherwise discard ~half of them).
ref_ids = Set.new
ref_count = 0
DuckDBUtil.query_each("SELECT value FROM run_metadata WHERE key = 'ref_count'"){ |l| ref_count = l.strip.to_i }
ref_ids.merge(node_ids.first(ref_count)) if ref_count > 0

## id2self from self_scores
id2self = {}
DuckDBUtil.query_each("SELECT seq_id, self_score FROM self_scores"){ |l|
	id, score = l.split("\t")
	id2self[id] = score.to_f
}

## summary_pre rows for kind='node'
id2out = Hash.new{ |h, i| h[i] = [] }
DuckDBUtil.query_each(
	"SELECT que_score_rank, sub_score_rank, que, sub, " +
	"que_score, sub_score, que_len, sub_len, " +
	"que_len_in_hit, sub_len_in_hit, " +
	"mean_idt_of_que, mean_idt_of_sub, " +
	"pct_que_len_in_hit, pct_sub_len_in_hit " +
	"FROM summary_pre WHERE kind = 'node' " +
	"ORDER BY node, que_score_rank"
){ |l|
	a = l.split("\t")[2..-1]  # drop ranks (match legacy "[2..-1]")
	que, sub = a[0..1]
	## Legacy filter: keep only [SHORTER => LONGER] direction, then write both sides.
	## In ref mode, input×ref HSPs exist in one direction only (ref is never queried);
	## skip the filter for those pairs so they aren't lost when input happens to be longer.
	## input×input pairs (both directions present) still get deduped.
	unless ref_ids.include?(que) || ref_ids.include?(sub)
		next if (id2order[que] <=> id2order[sub]) == 1
	end
	id2out[que] << a
	id2out[sub] << a if que != sub
}

## write combined TSV ready for COPY into summary_tsv
open(out_tsv, "w"){ |fout|
	node_ids.each{ |id|
		outs = []
		id2out[id].each{ |a|
			a = a.values_at(0..3, 8..11) # q, s, q_scr, s_scr, q_%idt, s_%idt, q_%len, s_%len
			query_is_id = a[0] == id
			targ  = query_is_id ? a[1] : a[0]
			p_idt = query_is_id ? a[4] : a[5]
			p_len = query_is_id ? a[6] : a[7]
			score = (a[2].to_i + a[3].to_i) * 0.5
			## SG normalises by self_score of the SHORTER sequence in the pair. In normal/2D modes
			## the shorter→longer filter ensures `a[0]` (que) is the shorter side, but ref mode
			## may keep input→ref HSPs where input is longer; pick the shorter id explicitly.
			shorter_id = ((id2order[a[0]] <=> id2order[a[1]]) <= 0) ? a[0] : a[1]
			s_scr = id2self[shorter_id]
			sg    = "%.4f" % [score / s_scr, 1].min            # max of SG should be 1
			outs << [targ, score, sg, p_idt, p_len]
		}
		outs.sort_by{ |a| -a[2].to_f }.each{ |a|
			fout.puts [id, id2info[a[0]], a[1..-1]].flatten*"\t"
		}
	}
}
