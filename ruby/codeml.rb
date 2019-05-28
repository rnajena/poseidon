#!/usr/local/bin/ruby

class Codeml

  attr_reader :mlcs

  CODEML_BIN = 'tools/paml4/codeml'

  CODON_FREQS = {:F61 => 0, :F1X4 => 1, :F3X4 => 2}
#  MODELS = [:M0, :M1a, :M2a, :M7, :M8]
  MODELS = [:M0, :M1a, :M2a, :M7, :M8, :M8a]

  def initialize(alignment, tree, dir, parameter_string)

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

      ctl_file = File.open("#{leave_path}/codeml.ctl",'w')
      codeml_ctls.push(write_ctl(ctl_file, alignment, tree, leave_path, codon_freq, model, ns_sites, fix_omega))
#      codeml_ctls.push(write_ctl_variation(ctl_file, alignment, tree, leave_path, codon_freq, model, ns_sites, fix_omega)) #TODO this was implemented to test different start omega values

    end

    codeml_mlcs = {:F61 => [], :F1X4 => [], :F3X4 => []}
    codeml_ctls.each do |ctl|
      dn = File.dirname(ctl)
      log_file = File.open("#{dn}/codeml.log",'w')
      freq = dn.split('/').reverse[1]
      codeml_mlcs[freq.intern].push("#{dn}/codeml.mlc")
      parameter_string << "codeml #{dn}/codeml.ctl\n"
      Process.fork do
        Dir.chdir(dn){
          log_file << `#{CODEML_BIN} codeml.ctl` unless File.exists?('codeml.mlc') #TODO
          log_file.close
        }
      end
    end

    puts 'I am waiting for the child processes running CODEML...'

    Process.waitall

    # combine all single codeml.mlc files to a final big codeml file
    @mlcs = combine_mlc(codeml_mlcs, dir)

  end

  ## for each frequency {F61, F1X4, F3X4} build one codeml mlc file
  def combine_mlc(codeml_mlcs, codeml_dir)
    mlcs = []
    codeml_mlcs.each do |freq, codeml_mlc_a|
      codeml_all_out = File.open("#{codeml_dir}/#{freq}/codeml.all.mlc",'w')
      codeml_mlc_a.each do |codeml_mlc|
        dn = File.dirname(codeml_mlc).split('/').reverse[0]
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
      mlcs.push(codeml_all_out.path)
    end
    mlcs
  end

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
    ctl_file << "seqfile  =  #{aln}                         * sequence data file name
treefile = #{tree}     * tree structure file name
outfile  = #{leave_path}/codeml.mlc                         * main result file name

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

  def write_ctl_variation(ctl_file, aln, tree, leave_path, codon_freq, model, ns_sites, fix_omega)
    ctl_file << "seqfile  =  #{aln}                         * sequence data file name
treefile = #{tree}     * tree structure file name
outfile  = #{leave_path}/codeml.mlc                         * main result file name

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
    omega = 0     * initial or fixed omega, for codons or codon-based AAs

       getSE = 0         * 0: don't want them, 1: want S.E.s of estimates
RateAncestor = 0         * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
  Small_Diff = .45e-6    * Default value.
   cleandata = 0         * remove sites with ambiguity data (1:yes, 0:no)
 fix_blength = 0         * 0: ignore, -1: random, 1: initial, 2: fixed\n"
    ctl_file.close
    ctl_file.path
  end

end
