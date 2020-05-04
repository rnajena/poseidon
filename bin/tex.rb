#!/usr/bin/env ruby

class Tex

  require 'bio'

  #PDFLATEX_BIN = 'pdflatex'
  #CHI2_BIN = '../tools/paml4/chi2'

  # check parameters for special symbols (that could fuck up LaTeX)
  SPECIAL_SYMBOLS = ['#', '$', '_', '%', 'ö', 'ä', 'ü', 'ß', '´', '&', '^', '{', '}', '\\', '"']

  attr_reader :output_file, :output_file_gapped, :gap_start2gap_length, :significance_m8, :significance_m2a, :significance_m8a

  def initialize(input, output, codon_freq, fragment_pos, title, reference_species, reference_aln, internal2input_species, aln_aa, frag_is_significant, chi2_bin)

    @chi2_bin = chi2_bin

    @internal2input_species = internal2input_species

    puts "\n%%%%%%%%%%%%%%%\n%% #{input}"

    # read in first (or reference) species aa sequence
    @reference_seq = nil
    @reference_id = nil
    save_reference_seq = nil; save_reference_id = nil
    reference_aln_length = 0
    Bio::FastaFormat.open(reference_aln).each do |e|
      id = e.definition.chomp
      reference_aln_length = e.seq.length
      if id == reference_species
        @reference_seq = e.seq.chomp
        @reference_id = id
        break
      end
      # just read the first species
      save_reference_seq = e.seq.chomp unless save_reference_seq
      save_reference_id = e.definition.chomp unless save_reference_id
    end
    unless @reference_id
      @reference_seq = save_reference_seq
      @reference_id = save_reference_id
    end

    frag_start = nil; frag_stop = nil
    if fragment_pos
      frag_start = fragment_pos.split('-')[0].to_i
      frag_stop = fragment_pos.split('-')[1].to_i
    end
    puts "Adjust positions because is a fragment: #{frag_start}-#{frag_stop}" if fragment_pos

    number_of_species = nil
    length_of_aln = nil

    read_m8_summary = false
    read_m8_beb = false

    percentage_sites_m8 = nil
    average_omega_of_sites_m8 = nil

    read_m2_summary = false
    read_m2_beb = false
    percentage_sites_m2 = nil
    average_omega_of_sites_m2 = nil

    beb_m8_entries = []
    beb_m2_entries = []

    input_file = File.open(input,'r')
    input_file.each do |line|

      number_of_species = line.split(' ')[0] if line.start_with?(' ') && !number_of_species
      length_of_aln = line.split(' ')[1] if line.start_with?(' ') && !length_of_aln

      read_m8_summary = true if line.start_with?('Model 8: beta')
      read_m2_summary = true if line.start_with?('Model 2: PositiveSelection')

      percentage_sites_m8 = line.split(' ') if line.start_with?('p:') && read_m8_summary && !percentage_sites_m8
      average_omega_of_sites_m8 = line.split(' ') if line.start_with?('w:') && read_m8_summary && !average_omega_of_sites_m8

      percentage_sites_m2 = line.split(' ') if line.start_with?('p:') && read_m2_summary && !percentage_sites_m2
      average_omega_of_sites_m2 = line.split(' ') if line.start_with?('w:') && read_m2_summary && !average_omega_of_sites_m2

      read_m8_beb = true if read_m8_summary && line.start_with?('Bayes Empirical Bayes (BEB) analysis')
      read_m2_beb = true if read_m2_summary && line.start_with?('Bayes Empirical Bayes (BEB) analysis')

      if read_m8_beb && line.start_with?(' ')
        s = line.split(' ')
        if s.length == 6 && s[0].to_i != 0
          beb_m8_entries.push(BebEntry.new(s, frag_start))
        end
      end

      if read_m2_beb && line.start_with?(' ')
        s = line.split(' ')
        if s.length == 6 && s[0].to_i != 0
          beb_m2_entries.push(BebEntry.new(s, frag_start))
        end
      end

    end
    input_file.close
    puts "read in #{beb_m8_entries.size} BEB entries for M8."

    m8_percent = (percentage_sites_m8[percentage_sites_m8.length-1].to_f * 100).round(2)
    m8_average_omega = (average_omega_of_sites_m8[average_omega_of_sites_m8.length-1].to_f).round(2)
    #puts "% sites with omega>1 -->  #{m8_percent}%, omega (dN/dS) = #{m8_average_omega}"

    m2_percent = (percentage_sites_m2[percentage_sites_m2.length-1].to_f * 100).round(2)
    m2_average_omega = (average_omega_of_sites_m2[average_omega_of_sites_m2.length-1].to_f).round(2)


    ## get log likelihoods to calculate significane values for M7 vs M8 and M1a vs M2a
    dof_m7 = nil; likelihood_m7 = nil; dof_m8 = nil; likelihood_m8 = nil
    dof_m8a = nil; likelihood_m8a = nil
    dof_m1a = nil; likelihood_m1a = nil; dof_m2a = nil; likelihood_m2a = nil
    model_count = 0
    `grep lnL #{input}`.split("\n").each do |lnl_line|
      lnl_entries = lnl_line.split(' ')
      dof = lnl_entries[3].split(')')[0].to_i
      likelihood = lnl_entries[4].to_f
      model_count += 1

      case model_count
        when 2
          dof_m1a = dof; likelihood_m1a = likelihood
        when 3
          dof_m2a = dof; likelihood_m2a = likelihood
        when 4
          dof_m7 = dof; likelihood_m7 = likelihood
        when 5
          dof_m8 = dof; likelihood_m8 = likelihood
        when 6
          dof_m8a = dof; likelihood_m8a = likelihood
      end
    end

    #puts "Model M1a: d.o.f. = #{dof_m1a}; log likelihood = #{likelihood_m1a}"
    #puts "Model M2a: d.o.f. = #{dof_m2a}; log likelihood = #{likelihood_m2a}"
    #puts "Model M7: d.o.f. = #{dof_m7}; log likelihood = #{likelihood_m7}"
    #puts "Model M8: d.o.f. = #{dof_m8}; log likelihood = #{likelihood_m8}"
    #puts "Model M8a: d.o.f. = #{dof_m8a}; log likelihood = #{likelihood_m8a}"

    # calculate LRT
    lrt_m7_m8 = (2* (likelihood_m8 - (likelihood_m7))).abs
    lrt_m1a_m2a = (2* (likelihood_m2a - (likelihood_m1a))).abs
    lrt_m8a_m8 = (2* (likelihood_m8 - (likelihood_m8a))).abs
    #puts "LRT: #{lrt}"
    #  compute difference in number of free parameters
    dof_m7_m8 = dof_m8 - dof_m7
    dof_m1a_m2a = dof_m2a - dof_m1a
    dof_m8a_m8 = dof_m8 - dof_m8a
    #puts "Degrees of freedom: #{dof}"
    # chi2 p --> #{dof} #{lrt}

    pvalue_m7_m8 = chi2(lrt_m7_m8, dof_m7_m8)
    pvalue_m1a_m2a = chi2(lrt_m1a_m2a, dof_m1a_m2a)
    pvalue_m8a_m8 = chi2(lrt_m8a_m8, dof_m8a_m8)

    @significance_m8 = LRT.new(lrt_m7_m8, pvalue_m7_m8, likelihood_m7, likelihood_m8, m8_percent, m8_average_omega)
    @significance_m2a = LRT.new(lrt_m1a_m2a, pvalue_m1a_m2a, likelihood_m1a, likelihood_m2a, m2_percent, m2_average_omega)
    @significance_m8a = LRT.new(lrt_m8a_m8, pvalue_m8a_m8, likelihood_m8a, likelihood_m8, m8_percent, m8_average_omega)

    title_save = ''
    title.scan(/./).each do |title_char|
      if SPECIAL_SYMBOLS.include?(title_char)
        title_save << ' '
      else
        title_save << title_char
      end
    end

    # build tex for M7 vs M8
    tex_m7_m8 = "\\documentclass[10pt,a4paper,oneside]{article}\n\\usepackage{multirow}\n\\usepackage{booktabs}\n\\usepackage{longtable}\n\\setlength\\LTcapwidth{1.4\\textwidth}\n\\setlength\\LTleft{0pt}\n\\setlength\\LTright{0pt}\n\n\\usepackage[top=1in, bottom=1.25in, left=1.0in, right=1.25in]{geometry}\n\n\\begin{document}\n\n"
    tex_m7_m8 << build_tex_table(beb_m8_entries, lrt_m7_m8, pvalue_m7_m8, number_of_species, length_of_aln, m8_percent, m8_average_omega, frag_start, frag_stop, codon_freq, title_save, nil, frag_is_significant, 'M7', 'M8')
    tex_m7_m8 << "\n\\end{document}\n"

    # build tex for M8a vs M8
    tex_m8a_m8 = "\\documentclass[10pt,a4paper,oneside]{article}\n\\usepackage{multirow}\n\\usepackage{booktabs}\n\\usepackage{longtable}\n\\setlength\\LTcapwidth{1.4\\textwidth}\n\\setlength\\LTleft{0pt}\n\\setlength\\LTright{0pt}\n\n\\usepackage[top=1in, bottom=1.25in, left=1.0in, right=1.25in]{geometry}\n\n\\begin{document}\n\n"
    tex_m8a_m8 << build_tex_table(beb_m8_entries, lrt_m8a_m8, pvalue_m8a_m8, number_of_species, length_of_aln, m8_percent, m8_average_omega, frag_start, frag_stop, codon_freq, title_save, nil, frag_is_significant, 'M8a', 'M8')
    tex_m8a_m8 << "\n\\end{document}\n"

    # build tex for M1a vs M2a
    tex_m1a_m2a = "\\documentclass[10pt,a4paper,oneside]{article}\n\\usepackage{multirow}\n\\usepackage{booktabs}\n\\usepackage{longtable}\n\\setlength\\LTcapwidth{1.4\\textwidth}\n\\setlength\\LTleft{0pt}\n\\setlength\\LTright{0pt}\n\n\\usepackage[top=1in, bottom=1.25in, left=1.0in, right=1.25in]{geometry}\n\n\\begin{document}\n\n"
    tex_m1a_m2a << build_tex_table(beb_m2_entries, lrt_m1a_m2a, pvalue_m1a_m2a, number_of_species, length_of_aln, m2_percent, m2_average_omega, frag_start, frag_stop, codon_freq, title_save, nil, frag_is_significant, 'M1a', 'M2a')
    tex_m1a_m2a << "\n\\end{document}\n"

    @output_file = File.open(output, 'w')
    @output_file << tex_m7_m8
    @output_file.close

    output_file_m8a_m8 = File.open(output.sub('.tex','.M8a_vs_M8.tex'), 'w')
    output_file_m8a_m8 << tex_m8a_m8
    output_file_m8a_m8.close

    output_file_m1a_m2a = File.open(output.sub('.tex','.M1a_vs_M2a.tex'), 'w')
    output_file_m1a_m2a << tex_m1a_m2a
    output_file_m1a_m2a.close

    # PDFTEX
    if (42 == 0)
    [@output_file, output_file_m8a_m8, output_file_m1a_m2a].each do |file|
      Process.fork do
        Dir.chdir(File.dirname(file)){
          `#{PDFLATEX_BIN} #{file.path}`
          `#{PDFLATEX_BIN} #{file.path}`
        }
      end
      Process.waitall
    end
    end

    ################################
    ## BUILD GAPPED OUTPUT

    # get gap positions from amino acid alignment
    gap_indices = []
    Bio::FastaFormat.open(aln_aa).each do |entry|
      arr = entry.seq.scan(/./)
      gap_indices += arr.each_index.select{|i| arr[i] == '-'}
    end
    gap_indices.sort!.uniq!

    gap_start2gap_length = {}

    gap_start = -1
    gap_end = -1
    i = 0
    gap_indices.each do |gap_pos|

      gap_start = gap_pos if gap_start == -1
      gap_end = gap_pos

      if gap_start != -1 && (gap_pos+1) != gap_indices[i+1]
        # we found the end of a gap, save the gap and reset
        gap_start2gap_length[gap_start+1] = (gap_end - gap_start + 1)
        gap_start = -1
        gap_end = -1
      end

      i += 1
    end
    if gap_start != -1
      gap_start2gap_length[gap_start+1] = (gap_end - gap_start + 1)
    end

    @gap_start2gap_length = gap_start2gap_length

    tex_gapped_m7_m8 = "\\documentclass[10pt,a4paper,oneside]{article}\n\\usepackage{multirow}\n\\usepackage{booktabs}\n\\usepackage{longtable}\n\\setlength\\LTcapwidth{1.4\\textwidth}\n\\setlength\\LTleft{0pt}\n\\setlength\\LTright{0pt}\n\n\\usepackage[top=1in, bottom=1.25in, left=1.0in, right=1.25in]{geometry}\n\n\\begin{document}\n\n"
    tex_gapped_m7_m8 << build_tex_table(beb_m8_entries, lrt_m7_m8, pvalue_m7_m8, number_of_species, length_of_aln, m8_percent, m8_average_omega, frag_start, frag_stop, codon_freq, title_save, gap_start2gap_length, frag_is_significant, 'M7', 'M8')
    tex_gapped_m7_m8 << "\n\\end{document}\n"

    tex_gapped_m8a_m8 = "\\documentclass[10pt,a4paper,oneside]{article}\n\\usepackage{multirow}\n\\usepackage{booktabs}\n\\usepackage{longtable}\n\\setlength\\LTcapwidth{1.4\\textwidth}\n\\setlength\\LTleft{0pt}\n\\setlength\\LTright{0pt}\n\n\\usepackage[top=1in, bottom=1.25in, left=1.0in, right=1.25in]{geometry}\n\n\\begin{document}\n\n"
    tex_gapped_m8a_m8 << build_tex_table(beb_m8_entries, lrt_m8a_m8, pvalue_m8a_m8, number_of_species, length_of_aln, m8_percent, m8_average_omega, frag_start, frag_stop, codon_freq, title_save, gap_start2gap_length, frag_is_significant, 'M8a', 'M8')
    tex_gapped_m8a_m8 << "\n\\end{document}\n"

    tex_gapped_m1a_m2a = "\\documentclass[10pt,a4paper,oneside]{article}\n\\usepackage{multirow}\n\\usepackage{booktabs}\n\\usepackage{longtable}\n\\setlength\\LTcapwidth{1.4\\textwidth}\n\\setlength\\LTleft{0pt}\n\\setlength\\LTright{0pt}\n\n\\usepackage[top=1in, bottom=1.25in, left=1.0in, right=1.25in]{geometry}\n\n\\begin{document}\n\n"
    tex_gapped_m1a_m2a << build_tex_table(beb_m2_entries, lrt_m1a_m2a, pvalue_m1a_m2a, number_of_species, length_of_aln, m2_percent, m2_average_omega, frag_start, frag_stop, codon_freq, title_save, gap_start2gap_length, frag_is_significant, 'M1a', 'M2a')
    tex_gapped_m1a_m2a << "\n\\end{document}\n"

    @output_file_gapped = File.open("#{output.sub('.tex','.gaps.tex')}",'w')
    @output_file_gapped << tex_gapped_m7_m8
    @output_file_gapped.close

    output_file_gapped_m8a_m8 = File.open("#{output.sub('.tex','.M8a_vs_M8.gaps.tex')}",'w')
    output_file_gapped_m8a_m8 << tex_gapped_m8a_m8
    output_file_gapped_m8a_m8.close

    output_file_gapped_m1a_m2a = File.open("#{output.sub('.tex','.M1a_vs_M2a.gaps.tex')}",'w')
    output_file_gapped_m1a_m2a << tex_gapped_m1a_m2a
    output_file_gapped_m1a_m2a.close

    if (42 == 0)
    [@output_file_gapped, output_file_gapped_m8a_m8, output_file_gapped_m1a_m2a].each do |file_gapped|
      Process.fork do
        Dir.chdir(File.dirname(file_gapped)){
          `#{PDFLATEX_BIN} #{file_gapped.path}`
          `#{PDFLATEX_BIN} #{file_gapped.path}`
        }
      end
      Process.waitall
    end
    end

  end

  def chi2(lrt, dof)

    p = `#{@chi2_bin} #{dof} #{lrt}`.split('=')[2].to_f.round(3)

    if p < 0.001
      '$<0.001$'
    else
      "$#{p}$"
    end
  end

  def build_tex_table(beb_entries, lrt, pvalue, num_of_species, length_of_aln, percent, average_omega, frag_start, frag_stop, codon_freq, title, gap_start2gap_length, frag_is_significant, m0, m1)

    pvalue_threshold = 0.1

    if pvalue.gsub('$','').to_f > pvalue_threshold
      percent = 'NA'
      average_omega = 'NA'
    end

    tex_table = "\\begin{longtable}{lcccccp{4cm}}\n"

    species_name = @reference_id
    if @internal2input_species[@reference_id]
      species_name = @internal2input_species[@reference_id]
    end
    species_name = species_name.gsub('_',' ')

    species_name_save = ''
    species_name.scan(/./).each do |char|
      if SPECIAL_SYMBOLS.include?(char)
        species_name_save << ' '
      else
        species_name_save << char
      end
    end


    if gap_start2gap_length
      tex_table << "\\caption{Results of the evolutionary analyses for positively selected sites for #{title}. P-values were achieved by performing chi-squared tests on twice the difference of the computed log likelihood values of the models disallowing (M7) or allowing (M8) $dN/dS>1$. The BEB column lists rapidly evolving sites with a $dN/dS>1$ and a posterior probability $>0.95$, determined by the Bayes Empirical Bayes implemented in \\texttt{Codeml}. Amino acids refer to \\emph{#{species_name_save}}. Note that INDELs and the stop codon were removed from the alignment prior to evolutionary analysis. Shown positions were mapped back to the \\textbf{alignment with gaps}.}\\\\\n"
    else
      tex_table << "\\caption{Results of the evolutionary analyses for positively selected sites for #{title}. P-values were achieved by performing chi-squared tests on twice the difference of the computed log likelihood values of the models disallowing (M7) or allowing (M8) $dN/dS>1$. The BEB column lists rapidly evolving sites with a $dN/dS>1$ and a posterior probability $>0.95$, determined by the Bayes Empirical Bayes implemented in \\texttt{Codeml}. Amino acids refer to \\emph{#{species_name_save}}. Note that INDELs and the stop codon were removed from the alignment prior to evolutionary analysis, so shown positions are based on the \\textbf{alignment without gaps}.}\\\\\n"
    end

    tex_table << "\\toprule\n"
    tex_table << "		\\bf Region & \\bf \\# species &  \\bf #{m0} vs #{m1} & \\bf #{m0} vs #{m1} & \\bf \\% sites with & \\bf avg($\\omega$) & \\bf #{m1} BEB \\\\
		 & & ($\\chi^2$) & p-value & $\\omega>1$ &  & ($PP>0.95/>\\textbf{0.99}$) \\\\\n\n\\endfirsthead\n\\multicolumn{7}{c}{{\\bfseries \\tablename\\ \\thetable{} -- continued from previous page}} \\\\ \n"
    tex_table << "\\toprule		\\bf Region & \\bf \\# species &  \\bf M7 vs M8 & \\bf M7 vs M8 & \\bf \\% sites with & \\bf avg($\\omega$) & \\bf M8 BEB \\\\
		 & & ($\\chi^2$) & p-value & $\\omega>1$ &  & ($PP>0.95/>\\textbf{0.99}$) \\\\\n\\midrule\n\\endhead\n\\multicolumn{7}{r}{{Continued on next page}} \\\\\n\\endfoot\n\\endlastfoot\n\n\\midrule\n"

    #tex_table = "\\begin{scriptsize}\n\\begin{table}\n\t\\begin{tabular}{lcccccp{4cm}}\n\t\t\\toprule\n\t\t\\bf Region & \\bf \\# Species & \\bf M7 vs M8 & \\bf M7 vs M8 & \\bf \\% sites with & \\bf avg($\\omega$) & \\bf M8 BEB \\\\ \n"
    #tex_table << "\t\t & & $\\chi^2$ & p-value & $\\omega>1$ &  & ($PP>0.95/>\\textbf{0.99}$) \\\\ \n\t\t\\midrule\n"

    tex_table << "\t\t \\multicolumn{7}{l}{\\textbf{#{codon_freq}}}  \\\\ \n\t\t \\midrule\n"

    if frag_start
      aln_start = frag_start
    else
      aln_start = 1
    end
    if frag_stop
      aln_stop = frag_stop
    else
      aln_stop = length_of_aln.to_i / 3
    end

    # adjust aln_start and aln_stop if we deal with a gapped output
    aln_start_adjustor = 0
    aln_stop_adjustor = 0
    if gap_start2gap_length
      gap_start2gap_length.each do |gap_start, gap_length|
        if gap_start < aln_start+aln_start_adjustor
          aln_start_adjustor += gap_length
        end
        if gap_start < aln_stop+aln_stop_adjustor
          aln_stop_adjustor += gap_length
        end
      end
    end

    #puts "#{m8_percent}\t#{m8_average_omega}"
    tex_table << "\t\taa #{aln_start+aln_start_adjustor}--#{aln_stop+aln_stop_adjustor} & #{num_of_species} & #{lrt.round(2)} & #{pvalue} & #{percent} & #{average_omega} & "

    bebs = []
    beb_entries.each do |beb_entry|
      # get the corresponding amino acid
      ref_aa = @reference_seq.scan(/./)[beb_entry.pos-1] # we have to do -1 because counting starts with 0 in the array!

      # adjust for the gapped alignment if necessary
      gap_pos_adjustor = 0
      current_pos = beb_entry.pos
      if gap_start2gap_length
        gap_start2gap_length.each do |gap_start, gap_length|
          if current_pos >= gap_start
            gap_pos_adjustor += gap_length
            current_pos += gap_length
          end
        end
      end

      if beb_entry.probability_omega > 0.99
        bebs.push("\\textbf{#{ref_aa}#{beb_entry.pos+gap_pos_adjustor}}")
      else
        if beb_entry.probability_omega > 0.95
          bebs.push("#{ref_aa}#{beb_entry.pos+gap_pos_adjustor}")
        end
      end
    end
    if pvalue.gsub('$','').to_f > pvalue_threshold
      tex_table << "NA \\\\ \n\t\t\\bottomrule \n"
    else
      if bebs.size > 0
        tex_table << "#{bebs.join('; ')} \\\\ \n\t\t\\bottomrule \n"
      else
        tex_table << "none \\\\ \n\t\t\\bottomrule \n"
      end
    end

    tex_table << "\\end{longtable}\n"

    tex_table
  end

end

class BebEntry

  attr_accessor :pos, :aminoacid, :probability_omega, :stars

  def initialize(values, frag_start)
    @pos = values[0].to_i
    @pos += frag_start - 1 if frag_start
    @aminoacid = values[1]
    @probability_omega = values[2].split('*')[0].to_f
    @stars = values[2].scan(/\*/).join('')

    #if frag_start == 184 && (@stars == '**' || @stars == '*')
    #  puts "#{values[0].to_i}\t#{@pos}\t#{@aminoacid}\t#{@probability_omega}\t#{@stars}"
    #end

  end

end

class LRT

  attr_accessor :lrt, :pvalue, :lnL0, :lnL1, :average_omega, :omega_percent

  def initialize(lrt, pvalue, lnL0, lnL1, omega_percent, average_omega)
    @lrt = lrt.to_f.round(2)
    @pvalue = pvalue
    @lnL0 = lnL0.to_f.round(2)
    @lnL1 = lnL1.to_f.round(2)
    @average_omega = average_omega
    @omega_percent = omega_percent
  end

end

#################################################
## INPUT 

mlc_file = ARGV[0]
freq = ARGV[1]
project_title = ARGV[2]
query_sequence_name = ARGV[3]
aln_aa_nogaps = ARGV[4]
internal2input_species_tsv = File.open(ARGV[5], 'r')
aln_aa = ARGV[6]
chi2_bin = ARGV[7]

fragment_pos = ARGV[8]
frag_is_significant = ARGV[9]

internal2input_species = {}
internal2input_species_tsv.each do |entry|
  s = entry.split("\t")
  internal2input_species[s[0]] = s[1].chomp
end

if fragment_pos
  Tex.new(mlc_file, mlc_file.sub('.mlc','.tex'), freq, fragment_pos, project_title, query_sequence_name, aln_aa_nogaps, internal2input_species, aln_aa, frag_is_significant, chi2_bin)
else
  Tex.new(mlc_file, mlc_file.sub('.mlc','.tex'), freq, nil, project_title, query_sequence_name, aln_aa_nogaps, internal2input_species, aln_aa, nil, chi2_bin)
end