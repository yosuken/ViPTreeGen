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
	puts "===== #{Time.now}"
	puts "===== step #{current} / #{total} (#{t.name}) -- #{status}"
	puts ""
	$stdout.flush
end
# }}} procedures


# {{{ default (run all tasks)
task :default do
	tasks = %w|
		01-1.prep_for_tblastx 01-2.tblastx            01-3.cat_and_rename_split_tblastx
		02-1.tblastx_filter   02-2.make_sbed_of_blast 02-3.make_summary_pre             02-4.make_self_tblastx 02-5.make_summary_tsv
		03-1.make_matrix      03-2.matrix_to_nj
	|
	NumStep = ENV["notree"] != "" ? tasks.size-1 : tasks.size
	tasks.each.with_index(1){ |task, idx|
		next if task =~ /^03-2/ and ENV["notree"] != "" # skip tree calculation if notree mode
		Rake::Task[task].invoke(idx)
	}
	# begin
	# rescue SystemExit => e
	# 	puts e.status
	# end
end
# }}} default (run all tasks)


# {{{ tasks 01
desc "01-1.prep_for_tblastx"
task "01-1.prep_for_tblastx", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	dir      = ENV["dir"]
	fin      = ENV["fin"]
	cutlen   = ENV["cutlen"].to_i # default: 100,000
	dbsize   = ENV["dbsize"]      # default: 200,000,000
	matrix   = ENV["matrix"]      # default: BLOSUM45
	evalue   = ENV["evalue"]      # default: 1e-2
	nthreads = "1"              
	jdir     = "#{dir}/batch/#{t.name.split(".")[0]}"; mkdir_p jdir
	odir     = "#{dir}/cat/all";    mkdir_p odir
	log      = "#{odir}/all.fasta.makeblastdb.log"
	fa       = "#{odir}/all.fasta"
	outs     = []
	max_target_seqs = 1_000_000

	## validate and make copy of input fasta
	open("#{odir}/all.len", "w"){ |flen|
		open(fa, "w"){ |fall|
			IO.read(fin).split(/^>/)[1..-1].each{ |ent|
				lab, *seq = ent.split("\n")
				seq  = seq.join
				len  = seq.size

				# take only label from comment line
				lab = lab.split[0]

				# label validation and modification
				lab = lab.gsub(/[^0-9a-zA-Z\_\-\.]/, "_")

				# [TODO] len  validation
				# raise if len  < XXX

				# make cat output
				flen.puts [lab, len]*"\t"
				fall.puts [">"+lab, seq.scan(/.{1,70}/)]

				# make node output
				# split fasta and make tblastx job
				n0dir = "#{dir}/node/#{lab}/seq";   mkdir_p n0dir
				n1dir = "#{dir}/node/#{lab}/blast"; mkdir_p n1dir

				_fsep = "#{n0dir}/#{lab}.fasta"
				open(_fsep, "w"){ |fsep|
					fsep.puts [">"+lab, seq.scan(/.{1,70}/)]
				}

				if len > cutlen
					n2dir = "#{n1dir}/split"; mkdir_p n2dir
					idx = 0
					while seq and seq.size > 0
						idx += 1
						subseq, seq = seq[0, cutlen], seq[cutlen..-1]
						_fspt = "#{n2dir}/#{idx}.fasta"
						open(_fspt, "w"){ |fspt|
							fspt.puts [">#{idx}", subseq.scan(/.{1,70}/)]
						}

						n3dir = "#{n2dir}/#{idx}"; mkdir_p n3dir
						_out  = "#{n3dir}/tblastx.out"
						_log  = "#{n3dir}/tblastx.log"
						outs << "tblastx -dbsize #{dbsize} -matrix #{matrix} -max_target_seqs #{max_target_seqs} -num_threads #{nthreads} \
						-evalue #{evalue} -outfmt 6 -db #{fa} -query #{_fspt} -out #{_out} 2>#{_log}".gsub(/\s+/, " ")
					end
				else
					_out = "#{n1dir}/tblastx.out"
					_log = "#{n1dir}/tblastx.log"
					outs << "tblastx -dbsize #{dbsize} -matrix #{matrix} -max_target_seqs #{max_target_seqs} -num_threads #{nthreads} \
					-evalue #{evalue} -outfmt 6 -db #{fa} -query #{_fsep} -out #{_out} 2>#{_log}".gsub(/\s+/, " ")
				end
			}
		}
	}

	## write batch file
	WriteBatch.call(outs, jdir, t)

	## makeblastdb
	sh "makeblastdb -dbtype nucl -hash_index -parse_seqids -in #{fa} -out #{fa} -title #{File.basename(fa)} 2>#{log}"
end
desc "01-2.tblastx"
task "01-2.tblastx", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	dir      = ENV["dir"]
	queue    = ENV["queue"]||""
	wtime    = ENV["wtime"]||"24:00:00"
	ncpus    = ENV["ncpus"]||""    
	nthreads = "1"
	mem      = nthreads.to_i * 12
	jdir     = "#{dir}/batch/01-1" # output from 01-1
	RunBatch.call(jdir, queue, nthreads, mem, wtime, ncpus)
end
desc "01-3.cat_and_rename_split_tblastx"
task "01-3.cat_and_rename_split_tblastx", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	dir      = ENV["dir"]
	cutlen   = ENV["cutlen"].to_i
	queue    = ENV["queue"]||"" 
	wtime    = ENV["wtime"]||"24:00:00"
	ncpus    = ENV["ncpus"]||""    
	nthreads = "1"
	mem      = nthreads.to_i * 12
	script   = "#{File.dirname(__FILE__)}/script/#{t.name}.rb"
	flen     = "#{dir}/cat/all/all.len"
	jdir     = "#{dir}/batch/#{t.name.split(".")[0]}"; mkdir_p jdir
	nodes    = IO.readlines(flen).inject({}){ |h, l| a=l.chomp.split("\t"); h[a[0]] = a[1].to_i; h }

	outs     = []
	nodes.each{ |node, len|
		n1dir = "#{dir}/node/#{node}/blast"
		if len > cutlen
			%w|tblastx|.each{ |program| 
				outs << "ruby #{script} #{node} #{program} #{flen} #{n1dir} #{cutlen}" 
			}
		end
	}
	WriteBatch.call(outs, jdir, t)
	RunBatch.call(jdir, queue, nthreads, mem, wtime, ncpus)
end
# }}} tasks 01


# {{{ tasks 02
desc "02-1.tblastx_filter"
task "02-1.tblastx_filter", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	dir      = ENV["dir"]
	cutlen   = ENV["cutlen"].to_i
	queue    = ENV["queue"]||""
	idt      = ENV["idt"]||"30"   # default: %idt > 30%
	aalen    = ENV["aalen"]||"30" # default: alignment length >= 30 aa.
	wtime    = ENV["wtime"]||"24:00:00"
	ncpus    = ENV["ncpus"]||""    
	nthreads = "1"
	mem      = nthreads.to_i * 12
	script   = "#{File.dirname(__FILE__)}/script/#{t.name}.rb"
	flen     = "#{dir}/cat/all/all.len"
	jdir     = "#{dir}/batch/#{t.name.split(".")[0]}"; mkdir_p jdir
	nodes    = IO.readlines(flen).inject({}){ |h, l| a=l.chomp.split("\t"); h[a[0]] = a[1].to_i; h }
	outs     = []

	nodes.each{ |node, len|
		n1dir = "#{dir}/node/#{node}/blast"
		if len > cutlen
			path = "#{n1dir}/split/*/tblastx.out2" # query names and postions are arranged by 01-3
		else
			path = "#{n1dir}/tblastx.out"
		end
		Dir[path].each{ |fin|
			fout = "#{fin.gsub(/2$/, "")}.filtered" # tblastx.out2 => tblastx.out.filtered
			outs << "ruby #{script} #{fin} #{fout} #{idt} #{aalen}"
		}
	}
	WriteBatch.call(outs, jdir, t)
	RunBatch.call(jdir, queue, nthreads, mem, wtime, ncpus)
end
desc "02-2.make_sbed_of_blast"
task "02-2.make_sbed_of_blast", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	dir      = ENV["dir"]
	cutlen   = ENV["cutlen"].to_i
	queue    = ENV["queue"]||""
	wtime    = ENV["wtime"]||"24:00:00"
	ncpus    = ENV["ncpus"]||""    
	nthreads = "1"
	mem      = nthreads.to_i * 12
	script   = "#{File.dirname(__FILE__)}/script/#{t.name}.rb"
	flen     = "#{dir}/cat/all/all.len"
	jdir     = "#{dir}/batch/#{t.name.split(".")[0]}"; mkdir_p jdir
	nodes    = IO.readlines(flen).inject({}){ |h, l| a=l.chomp.split("\t"); h[a[0]] = a[1].to_i; h }
	outs     = []

	nodes.each{ |node, len|
		n1dir = "#{dir}/node/#{node}/blast"
		%w|tblastx|.zip([".filtered"]){ |type, suffix|
			if len > cutlen
				path = "#{n1dir}/split/*/#{type}.out#{suffix}"
			else
				path = "#{n1dir}/#{type}.out#{suffix}"
			end
			Dir[path].each{ |fin|
				sbed = "#{File.dirname(fin)}/#{type}.sbed#{suffix}"
				#next if File.exist?(sbed) and !(File.zero?(sbed))
				outs << "ruby #{script} #{fin} #{sbed}" # --additional 12 --sort --scoring identity"
			}
		}
	}
	WriteBatch.call(outs, jdir, t)
	RunBatch.call(jdir, queue, nthreads, mem, wtime, ncpus)
end
desc "02-3.make_summary_pre"
task "02-3.make_summary_pre", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	dir      = ENV["dir"]
	cutlen   = ENV["cutlen"].to_i
	queue    = ENV["queue"]||""
	wtime    = ENV["wtime"]||"24:00:00"
	ncpus    = ENV["ncpus"]||""    
	nthreads = "1"
	mem      = nthreads.to_i * 12
	script   = "#{File.dirname(__FILE__)}/script/#{t.name}.rb"
	flen     = "#{dir}/cat/all/all.len"
	jdir     = "#{dir}/batch/#{t.name.split(".")[0]}"; mkdir_p jdir
	nodes    = IO.readlines(flen).inject({}){ |h, l| a=l.chomp.split("\t"); h[a[0]] = a[1].to_i; h }
	outs     = []

	nodes.each{ |node, len|
		n1dir = "#{dir}/node/#{node}/blast"
		%w|tblastx|.zip([".filtered"]){ |type, suffix|
			filename = "#{type}.sbed#{suffix}"
			fout  = "#{n1dir}/#{type}.summary.pre"
			outs << "ruby #{script} #{n1dir} #{fout} #{node} #{flen} #{filename}"
		}
	}
	# [TODO] rewrite the script without long computation
	WriteBatch.call(outs, jdir, t)
	RunBatch.call(jdir, queue, nthreads, mem, wtime, ncpus)
end
desc "02-4.make_self_tblastx"
task "02-4.make_self_tblastx", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	dir      = ENV["dir"]
	flen     = "#{dir}/cat/all/all.len"
	outs     = []
	bdir     = "#{dir}/cat/blast"; mkdir_p bdir

	%w|tblastx|.each{ |type|
		id2self = {}
		open("#{bdir}/#{type}.self.sum", "w"){ |fout|
			Dir["#{dir}/node/*/blast/#{type}.summary.pre"].each{ |fin|
				id = fin.split("/")[2]
				IO.readlines(fin)[1..-1].each{ |l|
					a = l.chomp.split("\t")
					id2self[id] = (a[4].to_i + a[5].to_i) * 0.5 if a[2] == a[3] # que == sub
				}
			}
			IO.readlines(flen).each{ |l|
				id = l.chomp.split("\t")[0]
				fout.puts [id, id2self[id]]*"\t"
			}
		}
	}
end
desc "02-5.make_summary_tsv"
task "02-5.make_summary_tsv", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	dir      = ENV["dir"]
	script   = "#{File.dirname(__FILE__)}/script/#{t.name}.rb"
	flen     = "#{dir}/cat/all/all.len"
	%w|tblastx|.each{ |type|
		sh "ruby #{script} #{flen} #{dir} #{type}" # output "#{dir}/node/*/blast/#{type}.summary.tsv"
	}
end
# }}} tasks 02


# {{{ tasks 03
desc "03-1.make_matrix"
task "03-1.make_matrix", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	dir     = ENV["dir"]
	flen    = "#{dir}/cat/all/all.len"
	rdir    = "#{dir}/result"; mkdir_p rdir

	ids  = []
	IO.readlines(flen).map{ |l| 
		ids << l.chomp.split("\t")[0]
	}

	similarities = Hash.new{ |h, i| h[i] = Hash.new(0.0) }
	%w|tblastx|.each{ |type|
		fout1   = open("#{rdir}/all.sim.matrix", "w")
		fout2   = open("#{rdir}/all.dist.matrix", "w")

		ids.each{ |id1|
			fin = "#{dir}/node/#{id1}/blast/#{type}.summary.tsv"
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
				out1 << ("%.4f" % similarities[id1][id2])
				out2 << ("%.4f" % (1 - similarities[id1][id2]))
				# [TODO] change this part in viptree
				# out << ( id1 == id2 ? 1 : "%.4f" % ( (similarities[id1][id2] + similarities[id2][id1]) * 0.5 ) )
			}
			fout1.puts out1*"\t"
			fout2.puts out2*"\t"
		}
		[fout1, fout2].each{ |fout| fout.close }
	}
end
desc "03-2.matrix_to_nj"
task "03-2.matrix_to_nj", ["step"] do |t, args|
	PrintStatus.call(args.step, NumStep, "START", t)
	dir      = ENV["dir"]
	notree   = ENV["notree"]
	method   = ENV["method"]||"bionj"
	rdir     = "#{dir}/result"
	flag     = "dist"
	script   = "#{File.dirname(__FILE__)}/script/#{t.name}.R"
	fin      = "#{rdir}/all.#{flag}.matrix"
	foutpref = "#{rdir}/all.#{method}"
	flog     = "#{rdir}/all.#{method}.makelog"
	sh "LANG=C Rscript --no-save --no-restore #{script} #{fin} #{foutpref} #{flag} #{method} >#{flog} 2>&1"
end
# }}} tasks 03
