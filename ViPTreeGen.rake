#
#  ViPTreeGen.rake - a command line tool for viral proteomic tree generation
#
#    Copyright: 2017 (C) Yosuke Nishimura (yosuke@kuicr.kyoto-u.ac.jp)
#    License: MIT license
#    Initial version: 2017-02-07
#


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

PrintStatus = lambda do |current, total, status, t|
	puts ""
	puts "\e[1;32m===== #{Time.now}\e[0m"
	puts "\e[1;32m===== step #{current} / #{total} (#{t.name}) -- #{status}\e[0m"
	puts ""
	$stdout.flush
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
					when "R"
						%|LANG=C R --version 2>&1|
					when "ape"
						%|LANG=C R --quiet --no-save --no-restore -e "packageVersion('ape')" 2>&1|
					when "phangorn"
						%|LANG=C R --quiet --no-save --no-restore -e "packageVersion('phangorn')" 2>&1|
					when "parallel"
						%{LANG=C parallel --version 2>&1 |head -n 1}
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
	tasks = %w|01-2.tblastx 01-3.cat_and_rename_split_tblastx 02-1.tblastx_filter 02-2.make_sbed_of_blast 02-3.make_summary_pre 02-4.make_self_tblastx|

	### add specific tasks
	if ENV["twoD"] == "" ## not 2D mode
		if ENV["notree"] == "" ## make tree
			tasks = %w|01-1.prep_for_tblastx| + tasks + %w|02-5.make_summary_tsv 03-1.make_matrix 03-2.matrix_to_nj|
		else ## only generate matrix
			tasks = %w|01-1.prep_for_tblastx| + tasks + %w|02-5.make_summary_tsv 03-1.make_matrix|
		end
	else ## 2D mode
		tasks = %w|01-1.2D.prep_for_tblastx| + tasks + %w|02-5.2D.make_summary_tsv 03-1.2D.make_matrix| ## add 2D tasks
	end

	### constants
	Odir     = ENV["dir"]
	Fin      = ENV["fin"]                  ## input file when 2D & normal mode
	Fin_q    = ENV["twoD"]                 ## query file only when 2D mode
	Flen     = "#{Odir}/cat/all/all.len"   ## create when 2D & normal mode
	Flen_q   = "#{Odir}/cat/all/query.len" ## create only when 2D mode
	Flen_i   = "#{Odir}/cat/all/input.len" ## create only when 2D mode
	Fa       = "#{Odir}/cat/all/all.fasta" ## create when 2D & normal mode

	Cutlen   = ENV["cutlen"].to_i ## default: 100,000
	DBsize   = ENV["dbsize"]      ## default: 200,000,000
	Matrix   = ENV["matrix"]      ## default: BLOSUM45
	Evalue   = ENV["evalue"]      ## default: 1e-2
	Nthreads = "1".to_i              
	Mem      = Nthreads * 12
	Qname    = ENV["queue"]||""
	Wtime    = ENV["wtime"]||"24:00:00"
	Ncpus    = ENV["ncpus"]||""    

	Max_target_seqs = 1_000_000

	TreeMethod = ENV["method"]||"bionj"

	### check version
	commands  = %w|tblastx makeblastdb R ape phangorn ruby|
	commands += %w|parallel| if Ncpus != ""
	CheckVersion.call(commands)

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
	PrintStatus.call(args.step, NumStep, "START", t)
	jdir     = "#{Odir}/batch/#{t.name.split(".")[0]}"; mkdir_p jdir
	odir     = "#{Odir}/cat/all"; mkdir_p odir
	log      = "#{odir}/all.fasta.makeblastdb.log"
	outs     = []
	puts ""

	## validate and make copy of input fasta
	fasta_str_q = IO.read(Fin_q)
	fasta_str_i = IO.read(Fin)

	# check if input is fasta format
	[fasta_str_q, fasta_str_i].zip(["query FASTA file", "input FASTA file"]){ |fasta_str, fasname|
		raise("\e[1;31mError:\e[0m #{fasname} might not be in FASTA format.") if fasta_str[0] != ">"
		raise("\e[1;31mError:\e[0m #{fasname} should include at least 3 seqeuences.") if fasname == "input FASTA file" and fasta_str.split(/^>/)[1..-1].size < 3
	}

	lab2seq_q = {}
	lab2seq_i = {}
	uniqlab   = {}
	[fasta_str_q, fasta_str_i].zip([lab2seq_q, lab2seq_i], ["query FASTA file", "input FASTA file"]){ |fasta_str, lab2seq, fasname|
		fasta_str.split(/^>/)[1..-1].each{ |ent|
			lab, *seq = ent.split("\n")
			seq = seq.join.gsub(/\s/, "")
			len = seq.size

			# take only label from comment line
			_lab = lab.split(/\s+/)[0]

			# label validation and modification
			not_allowed_pattern = /[^a-zA-Z0-9\-\.\_]/
			lab = _lab.gsub(not_allowed_pattern, "_") # modify signs to underscores

			### change hyphen and dots in start/end to underscore (begins with ".", "-" --> linux system problem, end with "." --> tblastx subject name problem ["abc." -> "abc"]), end with "___" --> ape/phangorn bionj generation problem ["abc___" -> "abc_"]
			not_allowed_pattern_pre = /^[\.\-\_]+/
			not_allowed_pattern_suf = /[\.\-\_]+$/
			not_allowed_pattern_dot = /(\.)[\.]+/
			not_allowed_pattern_hyp = /(\-)[\-]+/
			not_allowed_pattern_und = /(\_)[\_]+/
			## remove "." to "_" if the position is in the start / end  (issue in tblastx and tree generation step)
			## remove sequential "." or "-" or "_" (issue in bionj generation)
			lab = lab.gsub(not_allowed_pattern_pre, "").gsub(not_allowed_pattern_suf, "").gsub(not_allowed_pattern_dot, '\1').gsub(not_allowed_pattern_hyp, '\1').gsub(not_allowed_pattern_und, '\1')

			### label modification log
			puts "### [!] #{_lab} is changed to #{lab}" if _lab != lab

			### detect sequence format error
			raise("\e[1;31mError:\e[0m '#{lab}' is too short (length < 100 nt) included.") if len < 100
			# raise("\e[1;31mError:\e[0m not acceptable character is detected in sequence of '#{lab}' (allowed characters are 'ACGTRYKMSWBDHVN').") if seq =~ /[^ACGTRYKMSWBDHVN]/i
			raise("\e[1;31mError:\e[0m sequence name is not uniq. '#{lab}' is found multiple times in #{fasname}.") if lab2seq[lab]

			# store sequence name
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

			fque = "#{n0dir}/#{lab}.fasta"
			open(fque, "w"){ |_fque|
				_fque.puts [">"+lab, seq.scan(/.{1,70}/)]
			}

			if type == "input" ## make blastdb for each input sequence to compute self score
				sh "makeblastdb -dbtype nucl -in #{fque} -out #{fque} -title #{File.basename(fque)} 2>#{fque}.makeblastdb.log"
			end

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
					if type == "input"
						outs << "tblastx -dbsize #{DBsize} -matrix #{Matrix} -max_target_seqs #{Max_target_seqs} -num_threads #{Nthreads} \
						-evalue #{Evalue} -outfmt 6 -db #{fque} -query #{fspt} -out #{_out} 2>#{_log}".gsub(/\s+/, " ")
					else
						outs << "tblastx -dbsize #{DBsize} -matrix #{Matrix} -max_target_seqs #{Max_target_seqs} -num_threads #{Nthreads} \
						-evalue #{Evalue} -outfmt 6 -db #{Fa} -query #{fspt} -out #{_out} 2>#{_log}".gsub(/\s+/, " ")
					end
				end
			else
				_out = "#{n1dir}/tblastx.out"
				_log = "#{n1dir}/tblastx.log"
				if type == "input"
					outs << "tblastx -dbsize #{DBsize} -matrix #{Matrix} -max_target_seqs #{Max_target_seqs} -num_threads #{Nthreads} \
					-evalue #{Evalue} -outfmt 6 -db #{fque} -query #{fque} -out #{_out} 2>#{_log}".gsub(/\s+/, " ")
				else
					outs << "tblastx -dbsize #{DBsize} -matrix #{Matrix} -max_target_seqs #{Max_target_seqs} -num_threads #{Nthreads} \
					-evalue #{Evalue} -outfmt 6 -db #{Fa} -query #{fque} -out #{_out} 2>#{_log}".gsub(/\s+/, " ")
				end
			end
		}
	}

	## write batch file
	WriteBatch.call(outs, jdir, t)

	## makeblastdb
	sh "makeblastdb -dbtype nucl -in #{Fa} -out #{Fa} -title #{File.basename(Fa)} 2>#{log}"
end
desc "01-1.prep_for_tblastx"
task "01-1.prep_for_tblastx", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	jdir     = "#{Odir}/batch/#{t.name.split(".")[0]}"; mkdir_p jdir
	odir     = "#{Odir}/cat/all"; mkdir_p odir
	log      = "#{odir}/all.fasta.makeblastdb.log"
	fa       = "#{odir}/all.fasta"
	lab2seq  = {}
	outs     = []
	puts ""

	## validate and make copy of input fasta
	open(Flen, "w"){ |flen|
		open(fa, "w"){ |fall|
			fasta_str = IO.read(Fin)

			# check if input is fasta format
			raise("\e[1;31mError:\e[0m input FASTA file might not be in FASTA format.") if fasta_str[0] != ">"
			raise("\e[1;31mError:\e[0m input FASTA file should include at least 3 seqeuences.") if fasta_str.split(/^>/)[1..-1].size < 3

			uniqlab = {}
			fasta_str.split(/^>/)[1..-1].each{ |ent|
				lab, *seq = ent.split("\n")
				seq  = seq.join.gsub(/\s/, "")
				len  = seq.size

				# take only label from comment line
				_lab = lab.split(/\s+/)[0]

				# label validation and modification
				not_allowed_pattern = /[^a-zA-Z0-9\-\.\_]/
				lab = _lab.gsub(not_allowed_pattern, "_") # modify signs to underscores

				### change hyphen and dots in start/end to underscore (begins with ".", "-" --> linux system problem, end with "." --> tblastx subject name problem ["abc." -> "abc"]), end with "___" --> ape/phangorn bionj generation problem ["abc___" -> "abc_"]
				not_allowed_pattern_pre = /^[\.\-\_]+/
				not_allowed_pattern_suf = /[\.\-\_]+$/
				not_allowed_pattern_dot = /(\.)[\.]+/
				not_allowed_pattern_hyp = /(\-)[\-]+/
				not_allowed_pattern_und = /(\_)[\_]+/
				## remove "." to "_" if the position is in the start / end  (issue in tblastx and tree generation step)
				## remove sequential "." or "-" or "_" (issue in bionj generation)
				lab = lab.gsub(not_allowed_pattern_pre, "").gsub(not_allowed_pattern_suf, "").gsub(not_allowed_pattern_dot, '\1').gsub(not_allowed_pattern_hyp, '\1').gsub(not_allowed_pattern_und, '\1')

				### label modification log
				puts "### [!] #{_lab} is changed to #{lab}" if _lab != lab

				### detect sequence format error
				raise("\e[1;31mError:\e[0m '#{lab}' is too short (length < 100 nt) included.") if len < 100
				raise("\e[1;31mError:\e[0m sequence name is not uniq. '#{lab}' is found twice in #{fasname}.") if lab2seq[lab]
				### The rule below might be too strict --> do not use
				# raise("\e[1;31mError:\e[0m not acceptable character is detected in sequence of '#{lab}' (allowed characters are 'ACGTRYKMSWBDHVN').") if seq =~ /[^ACGTRYKMSWBDHVN]/i

				# store sequence name
				lab2seq[lab] = seq
				uniqlab[lab] = 1

				# make cat output
				flen.puts [lab, len]*"\t"
				fall.puts [">"+lab, seq.scan(/.{1,70}/)]

				# make node output
				# split fasta and make tblastx job
				n0dir = "#{Odir}/node/#{lab}/seq";   mkdir_p n0dir
				n1dir = "#{Odir}/node/#{lab}/blast"; mkdir_p n1dir

				fque = "#{n0dir}/#{lab}.fasta"
				open(fque, "w"){ |_fque|
					_fque.puts [">"+lab, seq.scan(/.{1,70}/)]
				}

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
						outs << "tblastx -dbsize #{DBsize} -matrix #{Matrix} -max_target_seqs #{Max_target_seqs} -num_threads #{Nthreads} \
						-evalue #{Evalue} -outfmt 6 -db #{Fa} -query #{fspt} -out #{_out} 2>#{_log}".gsub(/\s+/, " ")
					end
				else
					_out = "#{n1dir}/tblastx.out"
					_log = "#{n1dir}/tblastx.log"
					outs << "tblastx -dbsize #{DBsize} -matrix #{Matrix} -max_target_seqs #{Max_target_seqs} -num_threads #{Nthreads} \
					-evalue #{Evalue} -outfmt 6 -db #{Fa} -query #{fque} -out #{_out} 2>#{_log}".gsub(/\s+/, " ")
				end
			}
		}
	}

	## write batch file
	WriteBatch.call(outs, jdir, t)

	## makeblastdb
	sh "makeblastdb -dbtype nucl -hash_index -parse_seqids -in #{Fa} -out #{Fa} -title #{File.basename(Fa)} 2>#{log}"
end
desc "01-2.tblastx"
task "01-2.tblastx", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	jdir     = "#{Odir}/batch/01-1" # output from 01-1
	RunBatch.call(jdir, Qname, Nthreads, Mem, Wtime, Ncpus)
end
desc "01-3.cat_and_rename_split_tblastx"
task "01-3.cat_and_rename_split_tblastx", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	script   = "#{File.dirname(__FILE__)}/script/#{t.name}.rb"
	jdir     = "#{Odir}/batch/#{t.name.split(".")[0]}"; mkdir_p jdir
	outs     = []

	### set Nodes here
	if File.exist?(Flen_i) ## 2D mode
		Nodes    = IO.readlines(Flen_q).inject({}){ |h, l| a=l.chomp.split("\t"); h[a[0]] = a[1].to_i; h }
		Nodes_i  = IO.readlines(Flen_i).inject({}){ |h, l| a=l.chomp.split("\t"); h[a[0]] = a[1].to_i; h }
	else ## normal mode
		Nodes    = IO.readlines(Flen  ).inject({}){ |h, l| a=l.chomp.split("\t"); h[a[0]] = a[1].to_i; h }
		Nodes_i  = {}
	end

	[Nodes, Nodes_i].zip(%w|node input|){ |nodes, type|
		### Nodes_i --> for 2D mode only
		nodes.each{ |node, len|
			n1dir = "#{Odir}/#{type}/#{node}/blast"
			if len > Cutlen
				%w|tblastx|.each{ |program| 
					outs << "ruby #{script} #{node} #{program} #{Fa} #{n1dir} #{Cutlen}" 
				}
			end
		}
	}

	WriteBatch.call(outs, jdir, t)
	RunBatch.call(jdir, Qname, Nthreads, Mem, Wtime, Ncpus)
end
# }}} tasks 01


# {{{ tasks 02
desc "02-1.tblastx_filter"
task "02-1.tblastx_filter", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	idt      = ENV["idt"]||"30"   # default: %idt > 30%
	aalen    = ENV["aalen"]||"30" # default: alignment length >= 30 aa.
	script   = "#{File.dirname(__FILE__)}/script/#{t.name}.rb"
	jdir     = "#{Odir}/batch/#{t.name.split(".")[0]}"; mkdir_p jdir
	outs     = []

	[Nodes, Nodes_i].zip(%w|node input|){ |nodes, type|
		### Nodes_i --> for 2D mode only
		nodes.each{ |node, len|
			n1dir = "#{Odir}/#{type}/#{node}/blast"
			if len > Cutlen
				path = "#{n1dir}/split/*/tblastx.out2" # query names and postions are arranged by 01-3
			else
				path = "#{n1dir}/tblastx.out"
			end
			Dir[path].each{ |fin|
				fout = "#{fin.gsub(/2$/, "")}.filtered" # tblastx.out2 => tblastx.out.filtered
				outs << "ruby #{script} #{fin} #{fout} #{idt} #{aalen}"
			}
		}
	}
	WriteBatch.call(outs, jdir, t)
	RunBatch.call(jdir, Qname, Nthreads, Mem, Wtime, Ncpus)
end
desc "02-2.make_sbed_of_blast"
task "02-2.make_sbed_of_blast", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	script   = "#{File.dirname(__FILE__)}/script/#{t.name}.rb"
	jdir     = "#{Odir}/batch/#{t.name.split(".")[0]}"; mkdir_p jdir
	outs     = []

	[Nodes, Nodes_i].zip(%w|node input|){ |nodes, type|
		### Nodes_i --> for 2D mode only
		nodes.each{ |node, len|
			n1dir = "#{Odir}/#{type}/#{node}/blast"
			if len > Cutlen
				path = "#{n1dir}/split/*/tblastx.out.filtered"
			else
				path = "#{n1dir}/tblastx.out.filtered"
			end
			Dir[path].each{ |fin|
				sbed = "#{File.dirname(fin)}/tblastx.sbed.filtered"
				#next if File.exist?(sbed) and !(File.zero?(sbed))
				outs << "ruby #{script} #{fin} #{sbed}" # --additional 12 --sort --scoring identity"
			}
		}
	}
	WriteBatch.call(outs, jdir, t)
	RunBatch.call(jdir, Qname, Nthreads, Mem, Wtime, Ncpus)
end
desc "02-3.make_summary_pre"
task "02-3.make_summary_pre", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	script   = "#{File.dirname(__FILE__)}/script/#{t.name}.rb"
	jdir     = "#{Odir}/batch/#{t.name.split(".")[0]}"; mkdir_p jdir
	outs     = []

	[Nodes, Nodes_i].zip(%w|node input|){ |nodes, type|
		### Nodes_i --> for 2D mode only
		nodes.each{ |node, len|
			n1dir = "#{Odir}/#{type}/#{node}/blast"
			fout  = "#{n1dir}/tblastx.summary.pre"
			outs << "ruby #{script} #{n1dir} #{fout} #{node} #{Flen}"
		}
	}
	WriteBatch.call(outs, jdir, t)
	RunBatch.call(jdir, Qname, Nthreads, Mem, Wtime, Ncpus)
end
desc "02-4.make_self_tblastx"
task "02-4.make_self_tblastx", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	bdir     = "#{Odir}/cat/blast"; mkdir_p bdir
	outs     = []

	node2self = {}
	open("#{bdir}/tblastx.self.sum", "w"){ |fout|
		[Nodes, Nodes_i].zip(%w|node input|){ |nodes, type|
			### Nodes_i --> for 2D mode only
			nodes.each{ |node, len|
				fin = "#{Odir}/#{type}/#{node}/blast/tblastx.summary.pre"
				next unless File.exist?(fin)
				IO.readlines(fin)[1..-1].each{ |l|
					a = l.chomp.split("\t")
					node2self[node] = (a[4].to_i + a[5].to_i) * 0.5 if a[2] == a[3] # que == sub
				}
				fout.puts [node, node2self[node]]*"\t"
			}
		}
	}
end
desc "02-5.2D.make_summary_tsv" ### 2D mode
task "02-5.2D.make_summary_tsv", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	script   = "#{File.dirname(__FILE__)}/script/#{t.name}.rb"
	sh "ruby #{script} #{Flen} #{Odir}" # output "#{Odir}/node/*/blast/tblastx.summary.tsv"
end
desc "02-5.make_summary_tsv" ### normal mode
task "02-5.make_summary_tsv", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	script   = "#{File.dirname(__FILE__)}/script/#{t.name}.rb"
	sh "ruby #{script} #{Flen} #{Odir}" # output "#{Odir}/node/*/blast/tblastx.summary.tsv"
end
# }}} tasks 02


# {{{ tasks 03
desc "03-1.2D.make_matrix"
task "03-1.2D.make_matrix", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	rdir = "#{Odir}/result"; mkdir_p rdir
	qids = Nodes.keys
	tids = Nodes_i.keys

	similarities = Hash.new{ |h, i| h[i] = Hash.new(0.0) }
	qids.each{ |id1|
		fin = "#{Odir}/node/#{id1}/blast/tblastx.summary.tsv"
		IO.readlines(fin)[1..-1].each{ |l|
			sim, id2 = l.chomp.split("\t").values_at(3, 0)
			similarities[id1][id2] = sim.to_f
			similarities[id2][id1] = sim.to_f
		}
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
end
desc "03-1.make_matrix"
task "03-1.make_matrix", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	rdir = "#{Odir}/result"; mkdir_p rdir
	ids  = Nodes.keys

	similarities = Hash.new{ |h, i| h[i] = Hash.new(0.0) }
	fout1   = open("#{rdir}/all.sim.matrix", "w")
	fout2   = open("#{rdir}/all.dist.matrix", "w")

	ids.each{ |id1|
		fin = "#{Odir}/node/#{id1}/blast/tblastx.summary.tsv"
		IO.readlines(fin)[1..-1].each{ |l|
			sim, id2 = l.chomp.split("\t").values_at(3, 0)
			similarities[id1][id2] = sim.to_f
		}
	}

	[fout1, fout2].each{ |fout| fout.puts ["", ids]*"\t" }
	ids.each{ |id1|
		out1 = [id1]
		out2 = [id1]
		ids.each{ |id2|
			s = similarities[id1][id2]
			s = s == 0 ? "0" : (s == 1 ? "1" : "%.4f" % s)
			t = 1 - similarities[id1][id2]
			t = t == 0 ? "0" : (t == 1 ? "1" : "%.4f" % t)
			out1 << s
			out2 << t
		}
		fout1.puts out1*"\t"
		fout2.puts out2*"\t"
	}
	[fout1, fout2].each{ |fout| fout.close }
end
desc "03-2.matrix_to_nj"
task "03-2.matrix_to_nj", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	rdir     = "#{Odir}/result"
	flag     = "dist"
	script   = "#{File.dirname(__FILE__)}/script/#{t.name}.R"
	fin      = "#{rdir}/all.#{flag}.matrix"
	foutpref = "#{rdir}/all.#{TreeMethod}"
	flog     = "#{rdir}/all.#{TreeMethod}.makelog"
	sh "LANG=C Rscript --quiet --no-save --no-restore #{script} #{fin} #{foutpref} #{flag} #{TreeMethod} >#{flog} 2>&1"
end
# }}} tasks 03
