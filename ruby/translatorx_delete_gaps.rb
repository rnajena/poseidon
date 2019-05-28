#!/home/hoelzer/local/bin/ruby

require 'bio'

class TranslatorxDeleteGaps

  attr_reader :nt_aln_length

  TYPES = %w(nt aa)
  REVERSE_STOP_CODONS = %w(AAT GAT AGT AAU GAU AGU)

  def initialize(aln_dir, basename)
    TYPES.each do |type|
      gaps = nil
      seqs = {}

      nogaps = File.open("#{aln_dir}/#{basename}.#{type}_ali.nogaps.fasta",'w')
      remove(gaps, seqs, nogaps, aln_dir, basename, type)
    end
  end

  def remove(gaps, seqs, nogaps, aln_dir, bn, type)
    Bio::FastaFormat.open("#{aln_dir}/#{bn}.#{type}_ali.fasta").each do |entry|

      id = entry.definition
      seq = entry.seq

      # remove the last three nucleotieds or the last amino acid (stops)
      if type == 'aa'
        seq = seq.gsub('X','-').gsub('*','-')
      else
        # get last three chars that are not '-'
        if seq.reverse.start_with?('-')
          stop_codon = seq.reverse.scan(/-+.{3}/)[0].gsub('-','')
          seq = seq.reverse.sub(stop_codon,'---').reverse if REVERSE_STOP_CODONS.include?(stop_codon)
        else
          stop_codon = seq.reverse.scan(/.{3}/)[0]
          seq = seq.reverse.sub(stop_codon,'---').reverse if REVERSE_STOP_CODONS.include?(stop_codon)
        end
      end
  
      a = (0 ... seq.length).find_all { |i| seq[i,1] == '-' }
      if gaps
        gaps += a
      else
        gaps = a
      end
      seqs[id] = seq
    end
    gaps.uniq!

    seqs.each do |id, seq|
      nogaps << ">#{id}" << "\n"
      # remove gaps from seq
      seq_a = seq.scan(/.{0,1}/)
      nogap_seq = ''
      pos = 0
      seq_a.each do |nt|
        nogap_seq += nt unless gaps.include?(pos)
        pos += 1
      end
      nogaps << nogap_seq << "\n"
      @nt_aln_length = nogap_seq.length if type == 'nt'
    end
    nogaps.close
  end

end
