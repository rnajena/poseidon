#!/usr/bin/env ruby

# structure of an incoming lnc file: codeml_F3X4_M2a.mlc

MODELS = %w(M0 M1a M2a M7 M8 M8a)

codon_freq = ARGV[0]
codeml_mlcs = {}
codeml_mlcs[codon_freq] = []

## check if the mlc files are all empty and if so break
empty = false
Dir.glob("*.mlc").each do |mlc|
  empty = true if File.empty?(mlc)
end

if empty
  `touch dummy.all.mlc`
else
  ## we want to have a order of the models
  MODELS.each do |model|
    Dir.glob("*.mlc").each do |mlc|
      dn = mlc.split('_')[2].split('.')[0] # e.g. M2a
      if dn == model
        codeml_mlcs[codon_freq].push(mlc)
      end
    end
  end

## for each frequency {F61, F1X4, F3X4} build one codeml mlc file
#mlcs = []
codeml_mlcs.each do |freq, codeml_mlc_a|
  codeml_all_out = File.open("codeml_#{codon_freq}.all.mlc",'w')
  codeml_mlc_a.each do |codeml_mlc|
    puts codeml_mlc
    dn = codeml_mlc.split('_')[2].split('.')[0] # e.g. M2a
    mlc = File.open(codeml_mlc,'r')
    if dn == 'M0'
      mlc.each do |l|
        codeml_all_out << "\n\nModel 0: one-ratio\n\n" if l.start_with?('TREE')
        codeml_all_out << l
      end
      mlc.close
    else
      write = false
      mlc.each do |l|
        if l.start_with?('TREE')
          write = true
          case dn
            when 'M1a'
              codeml_all_out << "\n\nModel 1: NearlyNeutral (2 categories)\n\n"
            when 'M2a'
              codeml_all_out << "\n\nModel 2: PositiveSelection (3 categories)\n\n"
            when 'M7'
              codeml_all_out << "\n\nModel 7: beta (10 categories)\n\n"
            when 'M8'
              codeml_all_out << "\n\nModel 8: beta&w>1 (11 categories)\n\n"
            when 'M8a'
              codeml_all_out << "\n\nModel 8a: beta&w=1\n\n"
          end
        end
        codeml_all_out << l if write
      end
      mlc.close
    end
  end
  codeml_all_out.close
  #mlcs.push(codeml_all_out.path)
end
#mlcs

end