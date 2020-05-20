#!/usr/bin/env ruby

tex_summary_file = File.open(ARGV[0], 'w')
tex_summary_file_gaps = File.open(ARGV[1], 'w')

insignificant_frags = []
bp_tsv = File.open(ARGV[2],'r')
#stop_nt	stop_aa	significance	use_insignificance
#270	90	1	True
#549	183	1	True
#1948	650	1	True
frag_counter = 0
bp_tsv.each do |line|
  unless line.start_with?('#')
    s = line.split("\t")
    frag_counter += 1
    pvalue = s[2].to_f
    if pvalue > 0.1
      insignificant_frags.push(frag_counter)
    end
  end
end
bp_tsv.close

# parse all tex files into such a format:
# F61 => [codeml_F61.all.M8a_vs_M8.tex, codeml_F61.all.M1a_vs_M2a.tex], F3x4 => [...], ...
# and distinguish gaps and non-gap files
frequencies = %w(F61 F1X4 F3X4)

fragment_ids = []
Dir.glob("fragment_*_codeml_F3X4.all.M7_vs_M8.tex").each do |fragment_file|
  fragment_ids.push(fragment_file.split('_codeml_')[0])
end
fragment_ids.sort!

freq2tex_gapped = {}
frequencies.each do |freq|
  Dir.glob("codeml*.gaps.tex").each do |tex|
    # codeml_F61.all.gaps.tex
    next if tex.include?('M8a') || tex.include?('M2a') || !tex.include?(freq) # we only want to M7 vs M8 comparison
    codon_freq = File.basename(tex, '.gaps.tex').sub('codeml_','').split('.')[0]
    if codon_freq.start_with?('fragment')
      codon_freq = codon_freq.split('_')[2]
    end
    if freq2tex_gapped[codon_freq]
      freq2tex_gapped[codon_freq].push(tex) 
    else
      freq2tex_gapped[codon_freq] =  [tex] 
    end
  end
end
fragment_ids.each do |sorted_frag_id|
  frequencies.each do |freq|
    Dir.glob("*.gaps.tex").each do |tex|
      # codeml_F61.all.gaps.tex
      next if tex.include?('M8a') || tex.include?('M2a') || !tex.include?(freq) # we only want to M7 vs M8 comparison
      codon_freq = File.basename(tex, '.gaps.tex').sub('codeml_','').split('.')[0]
      if codon_freq.start_with?('fragment')
        codon_freq = codon_freq.split('_')[2]
      end
      if freq2tex_gapped[codon_freq]
        freq2tex_gapped[codon_freq].push(tex) if tex.start_with?(sorted_frag_id) || tex.start_with?('codeml_') && !freq2tex_gapped[codon_freq].include?(tex)
      else
        freq2tex_gapped[codon_freq] =  [tex] if tex.start_with?(sorted_frag_id) || tex.start_with?('codeml_')
      end
    end
  end
end

freq2tex = {}
frequencies.each do |freq|
  Dir.glob("codeml*.tex").each do |tex|
    next if tex.include?('.gaps.') || tex.include?('M8a') || tex.include?('M2a') || !tex.include?(freq) # we only want to M7 vs M8 comparison
    codon_freq = File.basename(tex, '.tex').sub('codeml_','').split('.')[0]
    if codon_freq.start_with?('fragment')
      codon_freq = codon_freq.split('_')[2]
    end
    if freq2tex[codon_freq]
      freq2tex[codon_freq].push(tex)
    else
      freq2tex[codon_freq] =  [tex]
    end
  end
end
fragment_ids.each do |sorted_frag_id|
  frequencies.each do |freq|
    Dir.glob("*.tex").each do |tex|
      # codeml_F61.all.gaps.tex
      # or
      # fragment_1_codeml_F61.all.M7_vs_M8.tex
      next if tex.include?('.gaps.') || tex.include?('M8a') || tex.include?('M2a') || !tex.include?(freq) # we only want to M7 vs M8 comparison
      codon_freq = File.basename(tex, '.tex').sub('codeml_','').split('.')[0]
      if codon_freq.start_with?('fragment')
        codon_freq = codon_freq.split('_')[2]
      end
      if freq2tex[codon_freq]
        freq2tex[codon_freq].push(tex) if tex.start_with?(sorted_frag_id) || tex.start_with?('codeml_') && !freq2tex[codon_freq].include?(tex)
      else
        freq2tex[codon_freq] =  [tex] if tex.start_with?(sorted_frag_id) || tex.start_with?('codeml_')
      end
    end
  end
end
puts freq2tex

def combine(freq2tex, tex_summary_file, insignificant_frags)
  caption = nil
  tex_header = ''
  tex_footer = ''
  tex_header << "\\documentclass[10pt,a4paper,oneside]{article}
  \\usepackage{multirow}
  \\usepackage{booktabs}
  \\usepackage{longtable}
  \\setlength\\LTcapwidth{1.3\\textwidth}
  \\setlength\\LTleft{0pt}
  \\setlength\\LTright{0pt}
  
  \\usepackage[top=1in, bottom=1.25in, left=1.0in, right=1.25in]{geometry}
  
  \\begin{document}
  
  \\begin{longtable}{lccccp{4cm}}\n"
  tex_footer << "\\toprule\n"
  tex_footer << "		\\bf Region &  \\bf M7 vs M8 & \\bf M7 vs M8 & \\bf \\% sites with & \\bf avg($\\omega$) & \\bf M8 BEB \\\\
    & ($\\chi^2$) & p-value & $\\omega>1$ &  & ($PP>0.95/>\\textbf{0.99}$) \\\\\n\n\\endfirsthead\n\\multicolumn{6}{c}{{\\bfseries \\tablename\\ \\thetable{} -- continued from previous page}} \\\\"
  tex_footer << "\\toprule		\\bf Region &  \\bf M7 vs M8 & \\bf M7 vs M8 & \\bf \\% sites with & \\bf avg($\\omega$) & \\bf M8 BEB \\\\
    & ($\\chi^2$) & p-value & $\\omega>1$ &  & ($PP>0.95/>\\textbf{0.99}$) \\\\\n\\midrule\n\\endhead\n\\multicolumn{6}{r}{{Continued on next page}} \\\\\n\\endfoot\n\\endlastfoot\n\n"
  freq2tex.each do |freq, tex_file_a|
    tex_footer << "		%% #{freq}\n"
    tex_footer << "	  \\midrule
    \\multicolumn{6}{l}{\\textbf{#{freq}}}  \\\\
    \\midrule\n"
    tex_file_a.each do |tex_file|
      frag = File.basename(tex_file).scan(/fragment_[0-9]+/)[0]
  
      if frag
        frag = frag.sub('ment_','')
        frag_count = frag.scan(/[0-9]+/)[0].to_i
        if insignificant_frags.include?(frag_count)
          frag = "#{frag}*"
        end
      else
        frag = 'full'
      end
      tex = File.open(tex_file,'r')
      tex.each do |tex_line|
        if tex_line.include?('aa') && !tex_line.include?('\caption')
          tex_split = tex_line.split('&')
          tex_split[1] = ''
          tex_split[0] = "#{frag} (#{tex_split[0].strip})"
          tex_tmp_split = []
          tex_split.each do |t|
            tex_tmp_split.push(t) unless t == ''
          end
          tex_string = tex_tmp_split.join(' & ')
          tex_footer << tex_string.gsub("\\","\\\\")
        end
        unless caption
          if tex_line.include?('\caption')
            caption = tex_line.gsub("\\","\\\\")
          end
        end
      end
      tex.close
    end
  end

  if insignificant_frags
    caption = caption.sub('}.}','}. Fragments arising from insignificant breakpoints (adjusted p-value $>$ 0.1) are marked with an asterisk.}') if insignificant_frags.size > 0
  end
  tex_summary_file << tex_header << caption.chomp << "\\\\\n" <<  tex_footer
  tex_summary_file << "\n\n"
  tex_summary_file << "\\bottomrule\n\n\\end{longtable}\n\\end{document}\n"
  tex_summary_file.close
  
end

combine(freq2tex, tex_summary_file, insignificant_frags)
combine(freq2tex_gapped, tex_summary_file_gaps, insignificant_frags)

