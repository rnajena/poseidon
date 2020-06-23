#!/usr/bin/env ruby

require 'bio'

    # check the final alignment for N's and if we have some, replace the full codon with ---
    # easiest way: get the X positions from the amino acid files and replace those positions
    
    #in
    nt_aln = Bio::FastaFormat.open("#{ARGV[0]}")
    aa_aln = Bio::FastaFormat.open("#{ARGV[1]}")
    html = File.open("#{ARGV[2]}",'r')

    #out
    nt_aln_tmp = File.open("#{ARGV[3]}",'w')
    aa_aln_tmp = File.open("#{ARGV[4]}",'w')
    html_tmp = File.open("#{ARGV[5]}",'w')

    log = File.open("#{ARGV[6]}",'r')

    # read in the log file
    logs = ''
    log.each do |line|
        logs << line
    end
    log.close
    `rm #{log.path}`
    log = File.open("#{ARGV[6]}",'w')
    log << logs

    @refactor = false
    aa_x_positions = {}
    aa_aln.each do |entry|
      seq = entry.seq
      id = entry.definition
      if seq.include?('X')
        next if seq.count('X') == 1 && seq.reverse.gsub('-','')[0] == 'X'
        @refactor = true
        #aa_x_positions[id] = seq.indices('X')
        aa_x_positions[id] = (0 ... seq.length).find_all { |i| seq[i,1] == 'X' }
        aa_x_positions[id].delete(seq.length-1) # because we dont want to remove the stop codon also encoded as X
        aa_x_positions.delete(id) if aa_x_positions[id].length == 0
      end
    end
    #aa_x_positions.uniq

    if aa_x_positions.values.join(' ').split(' ').uniq.length > 0
      puts "remove positions #{aa_x_positions} from the amino acid alignment."
      puts "\nPoSeiDon removed positions #{aa_x_positions.values.join(' ').split(' ').uniq} from the amino acid alignment because of ambiguous codons due to N bases in certain species.\n"
	    log << "\nPoSeiDon removed positions #{aa_x_positions.values.join(' ').split(' ').uniq} from the amino acid alignment because of ambiguous codons due to N bases in certain species.\n"
    end
    log.close

    if @refactor

      # refactor aa aln
      Bio::FastaFormat.open("#{ARGV[1]}").each do |entry|
        id = entry.definition
        seq_a = entry.seq.scan(/./)
        aa_aln_tmp << ">#{id}\n"
        i = 0
        seq_a.each do |aa|
          if aa_x_positions[id] && aa_x_positions[id].include?(i)
            aa_aln_tmp << '-'
          else
            aa_aln_tmp << aa
          end
          i += 1
        end
        aa_aln_tmp << "\n"
      end
      aa_aln_tmp.close

      # refactor nt aln
      nt_x_positions = {}
      aa_x_positions.each do |id, aa_positions|
        nt_x_positions[id] = []
        aa_positions.each do |aa_pos|
          nt_x_positions[id].push(aa_pos*3)
          nt_x_positions[id].push(aa_pos*3+1)
          nt_x_positions[id].push(aa_pos*3+2)
        end
      end
      nt_aln.each do |entry|
        id = entry.definition
        seq_a = entry.seq.scan(/./)
        nt_aln_tmp << ">#{id}\n"
        i = 0
        seq_a.each do |nt|
          if nt_x_positions[id] && nt_x_positions[id].include?(i)
            nt_aln_tmp << '-'
          else
            nt_aln_tmp << nt
          end
          i += 1
        end
        nt_aln_tmp << "\n"
      end
      nt_aln_tmp.close

      # refactor html
      check_lines = false
      i = 0
      species = ''
      html.each do |line|
        if line.include?('<table border=0><tr>')
          check_lines = true
        end

        if check_lines
          if line.include?('<tr><td>')
            # new species starts
            i = 0
            species = line.sub('</tr><tr><td>','').sub('</td>','').chomp
          end
          if line.include?('bgcolor')
            # check if we need to remove this codon line
            if aa_x_positions[species] && aa_x_positions[species].include?(i)
              html_tmp << "<td bgcolor=\"#FFFFFF\">---</td>\n"
            else
              html_tmp << line
            end
            i += 1
          else
            html_tmp << line
          end
        else
          html_tmp << line
        end
      end
      html_tmp.close
      html.close

      aa_aln.close
      nt_aln.close
    else

      nt_aln.each do |line|
        nt_aln_tmp << line
      end
      nt_aln.close; nt_aln_tmp.close

      aa_aln = Bio::FastaFormat.open("#{ARGV[1]}")
      aa_aln.each do |line|
        aa_aln_tmp << line
      end
      aa_aln.close; aa_aln_tmp.close
      
      html.each do |line|
        html_tmp << line
      end
      html.close; html_tmp.close

#      `cp #{aa_aln.path} #{aa_aln.path}.save`
#      `cp #{nt_aln.path} #{nt_aln.path}.save`
#      `cp #{html.path} #{html.path}.save`

#      `mv #{aa_aln_tmp.path} #{aa_aln.path}`
#      `mv #{nt_aln_tmp.path} #{nt_aln.path}`
#      `mv #{html_tmp.path} #{html.path}`

    end
