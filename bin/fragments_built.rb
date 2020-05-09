#!/usr/bin/env ruby

require 'bio'
require 'fileutils'

def build_fragments(aln, breakpoints)

    fragment_file_paths = []

    seqs = {}
    Bio::FastaFormat.open(aln).each do |entry|
      id = entry.definition
      seq = entry.seq
      seqs[id] = seq
    end
    puts "read in #{seqs.size} sequences and #{breakpoints.size} break points."

    seq_splits = {}
    seqs.keys.each do |id|
      seq_splits[id] = ''
    end
    break_counter = 0
    break_position_last = 0

    breakpoints.each do |breakpoint|
      break_counter += 1

      fragment_dir = "fragments/fragment_#{break_counter}/aln/"
      FileUtils.mkdir_p(fragment_dir)

      out = File.open("#{fragment_dir}/#{File.basename(aln)}",'w')
      fragment_file_paths.push(out.path)
      seqs.each do |id, seq|
        #out << ">#{id.chomp}  fragment=#{break_counter} breakpoint=#{(break_position_last.to_i+1).to_s}-#{breakpoint}\n"
        out << ">#{id.chomp}\n"
        seq_tmp = seq.scan(/.{1,#{breakpoint}}/)[0]
        seq_tmp = seq_tmp.sub(seq_splits[id],'')
        out << "#{seq_tmp.scan(/.{1,60}/).join("\n")}\n"
        seq_splits[id] += seq_tmp
      end
      break_position_last = breakpoint
      out.close
    end

    ## write out the rest
    #break_counter += 1

    #fragment_dir = "fragments/fragment_#{break_counter}/aln/"
    #FileUtils.mkdir_p(fragment_dir)

    #out = File.open("#{fragment_dir}/#{File.basename(aln)}",'w')
    #fragment_file_paths.push(out.path)
    #seqs.each do |id, seq|
    #  #out << ">#{id.chomp}  fragment=#{break_counter} breakpoint=#{(break_position_last.to_i+1).to_s}-#{seq.length}\n"
    #  out << ">#{id.chomp}\n"
    #  seq_tmp = seq
    #  seq_tmp = seq_tmp.sub(seq_splits[id],'')
    #  out << "#{seq_tmp.scan(/.{1,60}/).join("\n")}\n"
    #end
    #out.close

    # translate those nt fragments also to aa sequences
    fragment_file_paths.each do |fragment_nt_file|
      fragment_aa_file = File.open(fragment_nt_file.sub('.nt_ali','.aa_ali'),'w')
      in_file = Bio::FastaFormat.open(fragment_nt_file)
      in_file.each do |entry|
        id = entry.definition

        #seq = Bio::Sequence::auto(entry.seq)
        #aa_seq = seq.translate

        # create a DNA sequence
        seq = Bio::Sequence::NA.new(entry.seq.to_s)
        # translate to protein
        aa_seq = seq.translate
         
        fragment_aa_file << ">#{id.chomp}\n#{aa_seq.chomp}\n"
      end
      fragment_aa_file.close
      in_file.close
    end

    # write out the relative fragment file paths
    out = File.open('fragment_paths.txt','w')
    out << fragment_file_paths.join("\n") 
    out.close
    fragment_file_paths
end

############ Input

aln = ARGV[0]
breakpoints = [] # the nt positions of the bp
f = File.open(ARGV[1], 'r')
f.each do |line|
    pos = line.split("\t")[0].to_i
    breakpoints.push(pos) if pos > 0 
end
f.close

build_fragments(aln, breakpoints)
