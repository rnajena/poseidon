#!/usr/bin/env ruby

require 'fileutils'

dir = ARGV[0]
alignment = ARGV[1]
tree = ARGV[2]

CODON_FREQS = {:F61 => 0, :F1X4 => 1, :F3X4 => 2}
#  MODELS = [:M0, :M1a, :M2a, :M7, :M8]
MODELS = [:M0, :M1a, :M2a, :M7, :M8, :M8a]
NS_SITES_TO_MODEL = {0 => :M0, 1 => :M1a, 2 => :M2a, 3 => :M3, 7 => :M7, 8 => :M8}

def build_folder_structure(dir1)
  leave_folders = []
  Dir.mkdir(dir1) unless Dir.exists?(dir1)
  Dir.mkdir("#{dir1}/") unless Dir.exists?("#{dir1}/")
  CODON_FREQS.keys.each do |dir2|
    Dir.mkdir("#{dir1}/#{dir2}") unless Dir.exists?("#{dir1}/#{dir2}")
    MODELS.each do |model|
      Dir.mkdir("#{dir1}/#{dir2}/#{model}") unless Dir.exists?("#{dir1}/#{dir2}/#{model}")
      leave_folders.push("#{dir1}/#{dir2}/#{model}/")
    end
  end
  leave_folders
end


def write_ctl(ctl_file, aln, tree, leave_path, codon_freq, model, ns_sites, fix_omega)
  model_out = NS_SITES_TO_MODEL[ns_sites]
  model_out = 'M8a' if fix_omega == 1
  ctl_file << "seqfile  =  #{aln}                         * sequence data file name
treefile = #{tree}     * tree structure file name
outfile  = codeml_#{codon_freq}_#{model_out}.mlc                         * main result file name

noisy = 9       * 0,1,2,3,9: how much rubbish on the screen
verbose = 1       * 1: detailed output, 0: concise output
runmode = 0       * 0: user tree;  1: semi-automatic;  2: automatic 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

seqtype = 1     * 1:codons; 2:AAs; 3:codons-->AAs
CodonFreq = #{CODON_FREQS[codon_freq]}     * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
  clock = 0     * 0: no clock, unrooted tree, 1: clock, rooted tree
 aaDist = 0     * 0:equal, +:geometric; -:linear, {1-5:G1974,Miyata,c,p,v}
  model = #{model}     * models for codons:0:one, 1:b, 2:2 or more dN/dS ratios for branches.... or  0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85 5:T92, 6:TN93, 7:REV, 8:UNREST, 9:REVu; 10:UNRESTu

NSsites = #{ns_sites}       * 0:one w; 1:NearlyNeutral; 2:PositiveSelection; 3:discrete;4:freqs; 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;10:3normal
icode = 0       * 0:standard genetic code; 1:mammalian mt; 2-10:see below
Mgene = 0       * 0:rates, 1:separate; 2:pi, 3:kappa, 4:all

fix_kappa = 0     * 1: kappa fixed, 0: kappa to be estimated
  kappa = 2     * initial or fixed kappa
fix_omega = #{fix_omega}     * 1: omega or omega_1 fixed, 0: estimate
  omega = 1     * initial or fixed omega, for codons or codon-based AAs

     getSE = 0         * 0: don't want them, 1: want S.E.s of estimates
RateAncestor = 0         * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
Small_Diff = .45e-6    * Default value.
 cleandata = 0         * remove sites with ambiguity data (1:yes, 0:no)
fix_blength = 0         * 0: ignore, -1: random, 1: initial, 2: fixed\n"
  ctl_file.close
  ctl_file.path
end

## build the folder structure
leave_folders = build_folder_structure(dir)

# in each leave folder generate a specific ctl file to run codeml on
codeml_ctls = []
leave_folders.each do |leave_path|

  codon_freq = nil; model = 0; ns_sites = nil; fix_omega = 0

  leave_path.sub(dir, '').split('/').each do |parameter|
    if parameter.to_s.length > 1
          case parameter.intern
            when :F61
              codon_freq = :F61
            when :F1X4
              codon_freq = :F1X4
            when :F3X4
              codon_freq = :F3X4
            when :M0
              ns_sites = 0
            when :M1a
              ns_sites = 1
            when :M2a
              ns_sites = 2
            when :M3
              ns_sites = 3
            when :M7
              ns_sites = 7
            when :M8
              ns_sites = 8
            when :M8a
              ns_sites = 8
              fix_omega = 1
            else
              abort("I dont know this paramter to build config file from: #{parameter}")
          end
    end
  end

  model_out = NS_SITES_TO_MODEL[ns_sites]
  model_out = 'M8a' if fix_omega == 1
  ctl_file = File.open("#{leave_path}/codeml_#{codon_freq}_#{model_out}.ctl",'w')
  codeml_ctls.push(write_ctl(ctl_file, alignment, tree, leave_path, codon_freq, model, ns_sites, fix_omega))

end
