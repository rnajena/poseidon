#!/usr/bin/env ruby

require 'bio'

TYPES = %w(nt aa)
REVERSE_STOP_CODONS = %w(AAT GAT AGT AAU GAU AGU)

TYPES.each do |type|
  gaps = []
  seqs = {}

  if type == 'nt' 
    nogaps = File.open("#{ARGV[2]}",'w')
    aln = Bio::FastaFormat.open("#{ARGV[0]}")
  else
    nogaps = File.open("#{ARGV[3]}",'w')
    aln = Bio::FastaFormat.open("#{ARGV[1]}")
  end

  aln.each do |entry|

    id = entry.definition
    seq = entry.seq

    # remove the last three nucleotides or the last amino acid (stops)
    # the script before checked for internal stops already
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

    # check for gaps
    a = (0 ... seq.length).find_all { |i| seq[i,1] == '-' }
    a.each do |gap_pos|
      gaps.push(gap_pos) unless gaps.include?(gap_pos)
    end
    seqs[id] = seq
  end

  #gaps.uniq!
  aln.close

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