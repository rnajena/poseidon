#!/usr/bin/env ruby

require 'bio'

def build_aln_html_subfile(main_html, frag, bp_start, bp_end, main_aa_aln, gap_start2gap_length)

    puts "#{frag}\t#{bp_start}\t#{bp_end}"
    puts "GAPS: #{gap_start2gap_length}"

    sub_html = File.open("fragments/#{frag}/aln/#{File.basename(main_html)}",'w')
    sub_aa_aln = File.open("fragments/#{frag}/aln/#{File.basename(main_aa_aln)}",'w')

    main_nt_aln = main_aa_aln.sub('.aa_ali','.nt_ali')
    sub_nt_aln = File.open("fragments/#{frag}/aln/#{File.basename(main_aa_aln).sub('.aa_ali','.nt_ali')}",'w')

    write = true

    species_h = {}
    species_line = nil
    file = File.open(main_html,'r')
    file.each do |line|
      # here we reached the first species name
      if line.start_with?('</tr><tr><td>')
        write = false
        species_line = line
        species_h[line] = []
      end

      if write
        # write header and format stuff
        sub_html << line
      end

      if line.start_with?('<td bgcolor=') && species_line
        species_h[species_line].push(line)
      end
    end
    file.close

    # now loop over the species arrays and count the amount of gap columns between #{bp_start} and #{bp_end}
    gaps_until_start = 0
    gaps_until_end = 0

    gap_start2gap_length.each do |gap_start, gap_length|
      if (bp_start+gaps_until_start) > gap_start
        gaps_until_start += gap_length
      end
      if bp_end+gaps_until_end > gap_start
        gaps_until_end += gap_length
      end
    end

    puts "\t\tgaps until start: #{gaps_until_start}\tgaps until end: #{gaps_until_end}"

    #adjust break point positions according to gaps in the original alignment
    adj_bp_start = bp_start + gaps_until_start
    adj_bp_start -= 1 if bp_start != 0
    adj_bp_end = bp_end + gaps_until_end - 1

    species_h.each do |species_html_line, codon_array|
      sub_html << species_html_line
      (adj_bp_start..adj_bp_end).each do |adj_pos|
        sub_html << codon_array[adj_pos]
      end
    end

    sub_html << "</tr></table></div></div></body></html>\n"
    sub_html.close

    ## also adjust the aa_aln
    Bio::FastaFormat.open(main_aa_aln).each do |entry|
      id = entry.definition.chomp
      seq_a = entry.seq.chomp.scan(/./)
      sub_aa_aln << ">#{id}\n"
      (adj_bp_start..adj_bp_end).each do |adj_pos|
        sub_aa_aln << seq_a[adj_pos]
      end
      sub_aa_aln << "\n"
    end
    sub_aa_aln.close

    ## and the nt_aln
    Bio::FastaFormat.open(main_nt_aln).each do |entry|
      id = entry.definition.chomp
      seq_a = entry.seq.chomp.scan(/.{3}/)
      sub_nt_aln << ">#{id}\n"
      (adj_bp_start..adj_bp_end).each do |adj_pos|
        sub_nt_aln << seq_a[adj_pos]
      end
      sub_nt_aln << "\n"
    end
    sub_nt_aln.close

    aln_length_with_gaps = bp_end - bp_start + gaps_until_end

    #write out the gap-adjusted bp start and end for this fragment
    f = File.open("#{frag}_gap-adjusted_start_end.csv", 'w')
    f << "#{adj_bp_start},#{adj_bp_end}"
    f.close

    #[sub_html.path, sub_aa_aln.path, aln_length_with_gaps]
    puts "aln_length_with_gaps:#{aln_length_with_gaps}"
end


main_html = ARGV[0] # the original translatorx aa aln html

#stop_nt	stop_aa	significance	use_insignificance
#270	90	1	True
#549	183	1	True
#1948	650	1	True
frag = ARGV[1]
frag_count = 0
bp_start = 0
bp_end = 0
bp_tsv = File.open(ARGV[2],'r')
bp_tsv.each do |line|
    unless line.start_with?('#')
        s = line.split("\t")
        frag_count += 1
        bp = s[1].to_i
        if bp_end == 0
          bp_end = bp
        else
          bp_start = bp_end + 1
          bp_end = bp
        end
        break if "fragment_#{frag_count}" == frag
    end
end
bp_tsv.close

main_aa_aln = ARGV[3]

#gap_start,gap_length
#1,12
#21,1
#34,2
#43,1
#121,1
#125,1
#563,4
#580,1
gap_start2gap_length_csv = File.open(ARGV[4],'r')
gap_start2gap_length = {}
gap_start2gap_length_csv.each do |line|
    s = line.split(',')
    gap_start = s[0].to_i
    gap_length = s[1].to_i
    gap_start2gap_length[gap_start] = gap_length
end

frag_count = 0
f = File.open('aa_bp_with_gaps.csv' ,'w')
f << "#fragment,bp_aa_with_gaps,bp_aa_no_gaps,nr_aa_gaps\n"
bp_tsv = File.open(ARGV[2],'r')
bp_tsv.each do |line|
    unless line.start_with?('#')
        s = line.split("\t")
        frag_count += 1
        frag_gap_aa_length = 0
        aa_bp_no_gap = s[1].to_i
        gap_start2gap_length.each do |gap_start_aa_pos, gap_aa_length|
            puts "#{frag_gap_aa_length} += #{gap_aa_length} if #{aa_bp_no_gap} <= #{gap_start_aa_pos}"
            frag_gap_aa_length += gap_aa_length if aa_bp_no_gap > gap_start_aa_pos
        end
        f << "fragment_#{frag_count},#{frag_gap_aa_length+aa_bp_no_gap},#{aa_bp_no_gap},#{frag_gap_aa_length}\n"
    end
end
bp_tsv.close
f.close
##fragment,bp_aa_with_gaps
#fragment_1,106
#fragment_2,201
#fragment_3,673

build_aln_html_subfile(main_html, frag, bp_start, bp_end, main_aa_aln, gap_start2gap_length)