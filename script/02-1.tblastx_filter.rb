#
#  02-1.tblastx_filter.rb / ViPTreeGen - a command line tool for viral proteomic tree generation
#
#    Copyright: 2017 (C) Yosuke Nishimura (yosuke@kuicr.kyoto-u.ac.jp)
#    License: MIT license
#    Initial version: 2017-02-07
#
fin, fout, idt, aalen = ARGV
Idt = idt.to_i
AAlen = aalen.to_i

open(fout, "w"){ |fout|
	next unless File.exist?(fin)
	IO.readlines(fin).each{ |l|
		idt, aalen = l.chomp.split("\t").values_at(2, 3)
		next if idt.to_f  <= Idt
		next if aalen.to_i < AAlen
		fout.puts l
	}
}
