#!/home/hoelzer/local/bin/ruby

class Main

  require 'bio'
  require 'fileutils'

  require_relative('fasta_format_checker')
  require_relative('translatorx_delete_gaps')
  require_relative('phylo')
  require_relative('codeml')
  require_relative('html')
  require_relative('gard')
  require_relative('model_selection')
  require_relative('tex')
  require_relative('recombination_html')
  require_relative('codeml_html')
  require_relative('parameter_html')

  TRANSLATORX_BIN = 'tools/translatorx/translatorx.pl'
  POSEIDON_VERSION = '1.2'

  attr_reader :timestamp, :succeded_poseidon_run, :mail_notes, :is_recombination

  def initialize(dir, fasta, project_title, root_species, query_sequence_name, kh_insignificant_bp, user_mail, stamp)

    # if we need to refactor the alignment due to N's
    @refactor = false
    @is_recombination = false

    @succeded_poseidon_run = true

    if stamp
      @timestamp = stamp
    else
      @timestamp = Time.now.to_i
    end

    @threads = 2
    begin
      cpu = `lscpu -p | grep -v ^# | wc -l`.to_i
      load = `cut -d. -f 1 /proc/loadavg`.to_i
      free = cpu - load
      @threads = (free * 0.8).to_i
      @threads = 2 if @threads < 2
    rescue
      @threads = 2
    end

    @internal2input_species = {} # holds 'MYOTIS_DAVIDII' => 'Myotis_davidii'
    @input2internal_species = {} # holds 'Myotis_davidii => MYOTIS_DAVIDII'

    @nt_aln_length_nogaps = 0

    @root_species = []
    if root_species
      root_species.each do |species|
        @root_species.push(species.upcase.gsub('.','_').gsub('-','_').gsub(':','_').gsub(',','_').gsub(';','_'))
      end
    end
    @query_sequence_name = query_sequence_name.upcase.gsub('.','_').gsub('-','_').gsub(':','_').gsub(',','_').gsub(';','_')

    @dir = "#{dir}/results/#{@timestamp}"
    unless stamp
      # we have a generated time stamp, check if this time stamp already exist in some output folder
      while Dir.exists?("/mnt/fass2/poseidon-webserver-prost/results/#{@timestamp}") || Dir.exists?("/mnt/fass2/poseidon-webserver-mahlzeit/results/#{@timestamp}") || Dir.exists?("/mnt/fass2/poseidon-webserver-dessert/results/#{@timestamp}")
        @timestamp += 1
        @dir = "#{dir}/results/#{@timestamp}"
      end
    end
    Dir.mkdir(@dir) unless Dir.exists?(@dir)

    @parameters = {}

    @log = File.open("#{@dir}/#{@timestamp}.log",'w')
    @log << "###############################\n## START PoSeiDon: #{@timestamp}\n###############################\n\n"
    @parameters['full_aln'] = "###############################\n## START PoSeiDon: #{@timestamp}\n###############################\n\nPROJECT_DIR=/your/project/dir\n\n"

    puts "\nEstimated threads: #{@threads}\n"
    @log << "\nEstimated threads: #{@threads}\n"    
    
    bn = File.basename(fasta)
    bn = bn.split('.').reverse[1..bn.split('.').length].reverse.join('.')

    # tmp fasta file
    @fasta_file = File.open("#{@dir}/#{bn}.fna",'w')

    ## check fasta file and build a temporary file for further processing
    fasta_format_checker = FastaFormatChecker.new(fasta, @fasta_file, @internal2input_species, @input2internal_species, @log, @dir)
    fasta_format_checker.check(@query_sequence_name, user_mail, @timestamp, @root_species)
    @fasta_file.close
    @mail_notes = fasta_format_checker.mail_notes

    unless fasta_format_checker.continue
      @succeded_poseidon_run = false
      @log.close
      return
    end
    @parameters['full_aln'] << "#CALCULATIONS ON FULL ALIGNMENT...\n\n"

    ##############################################################
    ## 1) BUILD ALIGNMENT
    puts 'BUILD ALIGNMENT...'
    @log << "BUILD ALIGNMENT...\n"
    @parameters['full_aln'] << "#BUILD ALIGNMENT...\n"
    alignment_dir = "#{@dir}/aln"
    Dir.mkdir(alignment_dir) unless Dir.exists?(alignment_dir)
    alignment_out = "#{alignment_dir}/#{bn}"
    alignment(@fasta_file, alignment_out, @log)
    ##############################################################

    ##############################################################
    ## 2) REMOVE GAPS AND STOP CODONS FROM ALIGNMENT
    @log << "Remove gaps and stops from the alignment...\n"
    delete_gaps = TranslatorxDeleteGaps.new(alignment_dir, bn)
    @nt_aln_length_nogaps = delete_gaps.nt_aln_length

    aln_nt = "#{alignment_dir}/#{bn}.nt_ali.fasta"
    aln_aa = "#{alignment_dir}/#{bn}.aa_ali.fasta"
    aln_nt_nogaps = "#{alignment_dir}/#{bn}.nt_ali.nogaps.fasta"
    aln_aa_nogaps = "#{alignment_dir}/#{bn}.aa_ali.nogaps.fasta"
    aln_aa_html = "#{alignment_dir}/#{bn}.aa_based_codon_coloured.html"
    ##############################################################

    ##############################################################
    ## 3) PHYLOGENIES
    puts 'BUILD PHYLOGENY...'
    @log << "BUILD PHYLOGENY...\n"
    @parameters['full_aln'] << "#BUILD PHYLOGENY...\n"
    phylo_dir = "#{@dir}/phylo/"
    Dir.mkdir(phylo_dir) unless Dir.exists?(phylo_dir)

    out_dir = phylo_dir
    phylo_nt = Phylo.new(aln_nt_nogaps, nil, out_dir, :nt, @threads, @root_species, @log, @parameters['full_aln'])
    phylo_aa = Phylo.new(aln_nt_nogaps, nil, out_dir, :aa, @threads, @root_species, @log, @parameters['full_aln'])

    ##############################################################

    ##############################################################
    ## 4) AUTOMATIC MODEL SELECTION
    puts 'RUN MODEL SELECTION...'
    @log << "RUN MODEL SELECTION...\n"
    @parameters['full_aln'] << "#RUN MODELSELECTION...\n"
    model_out_dir = "#{@dir}/model_selection"
    FileUtils.mkdir_p(model_out_dir)
    model_selection = ModelSelection.new(aln_nt_nogaps, phylo_nt.tree_unshod, 4, 1, 0.05, model_out_dir, @threads, @parameters['full_aln'])
    nucleotide_bias_model = model_selection.model
    ##############################################################

    ##############################################################
    ## 5) GARD BREAKPOINT ANALYSIS
    puts 'RUN GARD...'
    @log << "RUN GARD...\n"
    @parameters['full_aln'] << "#RUN GARD...\n"
    gard_out_dir = "#{@dir}/gard/"
    Dir.mkdir(gard_out_dir) unless Dir.exists?(gard_out_dir)
    gard = Gard.new(aln_nt_nogaps, nucleotide_bias_model, '2', '3', "#{gard_out_dir}/gard", @threads, kh_insignificant_bp, @parameters['full_aln'])
    @mail_notes << gard.mail_notes if gard.mail_notes.length > 1
    puts 'GARD DONE.'

    breakpoints = gard.breakpoints
    gard_html_file = gard.gard_html_file

    ##############################################################

    ##############################################################
    ## 6) CODEML
    puts 'CODEML on full alignment...'
    @log << "CODEML on full alignment...\n"
    @parameters['full_aln'] << "#RUN CODEML ON FULL ALIGNMENT...\n"
    codeml_full_dir = "#{@dir}/codeml/"
    FileUtils.mkdir_p(codeml_full_dir)
    codeml = Codeml.new(aln_nt_nogaps, phylo_nt.tree_unshod, codeml_full_dir, @parameters['full_aln'])
    puts 'CODEML DONE on full alignment.'
    @log << "CODEML DONE on full alignment...\n"
    ## BUILD LATEX SUMMARY TABLE
    tex_full_aln = nil
    tex_objects = {'full_aln' => []}
    codeml.mlcs.each do |mlc_file|
      freq = File.dirname(mlc_file).split('/').reverse[0]
      tex_full_aln = Tex.new(mlc_file, mlc_file.sub('.mlc','.tex'), freq, nil, project_title, @query_sequence_name, aln_aa_nogaps, @internal2input_species, aln_aa, nil)
      tex_objects['full_aln'].push(tex_full_aln)
    end
    @log << "\tTex done on full alignment.\n"

    #add the full aln length to the breakpoints array to calculate also the last fragment, if we have fragments!
    aa_breakpoints = []
    breakpoints.keys.each do |nt_bp|
      aa_breakpoints.push(nt_bp/3)
    end
    aa_breakpoints.push(@nt_aln_length_nogaps/3 + 1) if breakpoints.keys.size > 0

    # 6.1) CODEML on FRAGMENTS, if there are any
    fragment_file_paths = nil
    phylo_fragments = nil
    if breakpoints.keys.size > 0

      @is_recombination = true

      @log << "Work on fragments...\n"
      @parameters['full_aln'] << "\n\n#BREAKPOINTS DETECTED, ALIGNMENT ACCORDINGLY SPLITTED\n"

      fragments_dir = "#{@dir}/fragments"
      Dir.mkdir(fragments_dir) unless Dir.exists?(fragments_dir)

      # build GARD fragments of alignment
      fragment_file_paths = build_fragments(aln_nt_nogaps, breakpoints.keys)
      sleep 5

      # for each fragment, build a new unrooted tree for CODEML
      phylo_fragments = {}
      fragment_file_paths.each do |fragment_file|
        frag = File.dirname(fragment_file).scan(/fragment_[0-9]+/)[0]
        @parameters[frag] = "#BREAKPOINTS DETECTED, ALIGNMENT ACCORDINGLY SPLITTED\n\n################################################\n#CALCULATIONS ON FRAGMENT #{frag.split('_')[1].chomp}\n\n#BUILD PHYLOGENY #{frag}\n"
        out_dir = "#{@dir}/fragments/#{frag}/phylo/"
        FileUtils.mkdir_p(out_dir)
        phylo_fragments[fragment_file] = (Phylo.new(fragment_file, nil, out_dir, :nt, @threads, @root_species, @log, @parameters[frag])) # [bats_frag1.fasta => PHYLO_OBJECT]
        phylo_fragments[fragment_file.sub('.nt_ali','.aa_ali')] = (Phylo.new(fragment_file.sub('.nt_ali','.aa_ali'), nil, out_dir, :aa, @threads, @root_species, @log, @parameters[frag])) # [bats_frag1.aa.fasta => PHYLO_OBJECT]
      end

      # check according to the initial breakpoint array if they are significant
      insignificant_frags = []

      bp_signi_counter = 0
      breakpoints.each do |bp, signi|
        bp_signi_counter += 1
        if signi == 1
          insignificant_frags.push(bp_signi_counter) unless insignificant_frags.include?(bp_signi_counter)
          insignificant_frags.push(bp_signi_counter+1) unless insignificant_frags.include?(bp_signi_counter+1)
        end
      end

      # now run CODEML for each fragment separately
      bp_pos = -1
      fragment_file_paths.each do |fragment_file|
        frag = File.dirname(fragment_file).scan(/fragment_[0-9]+/)[0]
        frag_counter = frag.scan(/[0-9]+/)[0].to_i
        codeml_dir = "#{@dir}/fragments/#{frag}/codeml/"
        tex_objects[frag] = []
        FileUtils.mkdir_p(codeml_dir)
        @parameters[frag] << "#RUN CODEML ON #{frag}\n"
        codeml = Codeml.new(fragment_file, phylo_fragments[fragment_file].tree_unshod, codeml_dir, @parameters[frag])
        ## BUILD LATEX SUMMARY TABLE for each fragment separately
        codeml.mlcs.each do |mlc_file|
          freq = File.dirname(mlc_file).split('/').reverse[0]
          frag_is_significant = true
          if frag == 'fragment_1'
            frag_is_significant = false if insignificant_frags.include?(frag_counter)
            tex = Tex.new(mlc_file, mlc_file.sub('.mlc','.tex'), freq, nil, "#{project_title} (Fragment #{frag.split('_')[1]})", @query_sequence_name, aln_aa_nogaps, @internal2input_species, aln_aa, frag_is_significant)
            tex_objects[frag].push(tex)
          else
            frag_is_significant = false if insignificant_frags.include?(frag_counter)
            bp_start = aa_breakpoints[bp_pos] + 1
            bp_end = aa_breakpoints[bp_pos+1]
            bp_end -= 1 if bp_end == aa_breakpoints.reverse[0]
            puts "#{bp_start}-#{bp_end}"
            tex = Tex.new(mlc_file, mlc_file.sub('.mlc','.tex'), freq, "#{bp_start}-#{bp_end}", "#{project_title} (Fragment #{frag.split('_')[1]})", @query_sequence_name, aln_aa_nogaps, @internal2input_species, aln_aa, frag_is_significant)
            tex_objects[frag].push(tex)
          end
        end
        bp_pos += 1
      end

    end

    frag_names = []
    if fragment_file_paths
      fragment_file_paths.each do |fragment_file|
        frag_names.push(File.dirname(fragment_file).scan(/fragment_[0-9]+/)[0])
      end
    end

    # refactor command paths
    puts @parameters.keys
    @parameters['full_aln'].gsub!(@dir,'$PROJECT_DIR')
    @parameters['full_aln'].gsub!('tools/','')
    @parameters['full_aln'].gsub!('translatorx/','')
    @parameters['full_aln'].gsub!('raxml/8.0.25/','')
    @parameters['full_aln'].gsub!('nw_utilities/','')
    @parameters['full_aln'].gsub!('//','/')
    frag_names.each do |frag|
      @parameters[frag].gsub!(@dir,'$PROJECT_DIR')
      @parameters[frag].gsub!('tools/','')
      @parameters[frag].gsub!('translatorx/','')
      @parameters[frag].gsub!('raxml/8.0.25/','')
      @parameters[frag].gsub!('nw_utilities/','')
      @parameters[frag].gsub!('//','/')
    end

    puts 'BUILD COMBINED TEX AND PDF FILE'
    @log << "BUILD COMBINED TEX AND PDF OUTPUT...\n"
    #combine all codeml TEX tables into one big table ('longtable'), publication-ready. Do this for the gapped and the ungapped alignment.
    #ungapped
    freq2tex = {:F61 => [], :F1X4 => [], :F3X4 => []}
    tex_objects.values.each do |object_a|
      object_a.each do |object|
        tex_file_path = object.output_file.path
        freq = File.dirname(tex_file_path).split('/').reverse[0]
        freq2tex[freq.intern].push(tex_file_path)
      end
    end
    tex_summary_file = File.open("#{@dir}/#{File.basename(@fasta_file,'.fna')}.tex",'w')
    build_combined_tex(freq2tex, tex_summary_file, insignificant_frags)

    #gapped
    freq2tex = {:F61 => [], :F1X4 => [], :F3X4 => []}
    tex_objects.values.each do |object_a|
      object_a.each do |object|
        tex_file_path = object.output_file.path.sub('.tex','.gaps.tex')
        freq = File.dirname(tex_file_path).split('/').reverse[0]
        freq2tex[freq.intern].push(tex_file_path)
      end
    end
    tex_summary_file_gapped = File.open("#{@dir}/#{File.basename(@fasta_file,'.fna')}.gaps.tex",'w')
    build_combined_tex(freq2tex, tex_summary_file_gapped, insignificant_frags)

    ##############################################################

    ##############################################################
    ## 7) HTML OUTPUT
    puts "BUILD HTML OUTPUT"
    @log << "BUILD HTML OUTPUT...\n"

    #breakpoint_pos_with_gaps = [] # holds the gap-adjusted breakpoint positions

    index_html_paths = {}
    index_html_paths['full_aln'] = '../full_aln/index.html'
    aa_breakpoints.size.times do |i|
      index_html_paths["fragment_#{i+1}"] = "../fragment_#{i+1}/index.html"
    end

    ## 7.0) build also HTML output for the fragments, if any
    colors = %w(#D2691E #0000CD #006400 #4B0082 #800000 #2E8B57 #00CED1 #808000 #FF69B4 #FF1493 #696969)

    adjusted_domain_pos = {}

    if fragment_file_paths

      bp_pos = 0
      bp_start = 0
      bp_end = 0

      aln_length_with_gaps = 0 # we need this to adjust the html first row column counting

      fragment_file_paths.each do |fragment_file|

        frag = File.dirname(fragment_file).scan(/fragment_[0-9]+/)[0]
        phylo_object = phylo_fragments[fragment_file]

        html_dir = "#{@dir}/html/"
        html_out = "#{@dir}/html/#{frag}/index.html"

        if @root_species
          tree_nt = phylo_object.tree_rooted
          tree_aa = phylo_object.tree_rooted.sub('.nt.','.aa.')
        else
          tree_nt = phylo_object.tree_corrected
          tree_aa = phylo_object.tree_corrected.sub('.nt.','.aa.')
        end
#        tree_nt = phylo_object.tree_corrected
#        tree_aa = phylo_object.tree_corrected

        #1) from the original translatorX aln_aa_html output, build a subfile with only the fragment sequence
        #2) aln_aa should be also a subfile
        #3) @fasta_file is a subfile (the fragment fasta)... should be fine... is the original file, so what.

        #1)
        bp = aa_breakpoints[bp_pos]
        if bp_end == 0
          bp_end = bp
        else
          bp_start = bp_end + 1
          bp_end = bp
        end
        results_a = build_aln_html_subfile(aln_aa_html, frag, bp_start, bp_end, aln_aa, tex_full_aln.gap_start2gap_length)
        aln_aa_html_sub = results_a[0]

        #2) here, we can just use the positions already calculated in 1) to adjust the aa_aln from translatorx
        aln_aa_sub = results_a[1]

        # colored bar for the fragment
        domain_frag_pos = {}
        fragment_count = frag.split('_')[1].to_i
        count_adjustor = nil
        if fragment_count > colors.size
          count_adjustor = colors.size
          domain_frag_pos["Fragment_#{fragment_count}".intern] = [[Range.new(bp_start,bp_end), colors[fragment_count-1-count_adjustor],true]]
        else
          domain_frag_pos["Fragment_#{fragment_count}".intern] = [[Range.new(bp_start,bp_end), colors[fragment_count-1],true]]
        end

        # adjust the domain positions to gaps we have in the alignment
        adjusted_domain_frag_pos = {}
        domain_frag_pos.each do |domain_label, domain_array|
          domain_array.each do |domain_entry|
            range = domain_entry[0]
            color = domain_entry[1]
            visible = domain_entry[2]

            bp_start = range.min
            bp_end = range.max

            gaps_until_start = 0
            gaps_until_end = 0

            tex_full_aln.gap_start2gap_length.each do |gap_start, gap_length|
              if (bp_start+gaps_until_start) > gap_start
                gaps_until_start += gap_length
              end
              if bp_end+gaps_until_end > gap_start
                gaps_until_end += gap_length
              end
            end
            bp_start = 1 if bp_start == 0
            adjusted_domain_frag_pos[domain_label] = [[Range.new(bp_start+gaps_until_start-1,bp_end+gaps_until_end-1), color, visible]]
            adjusted_domain_pos[domain_label] = [[Range.new(bp_start+gaps_until_start-1,bp_end+gaps_until_end-1), color, visible]] # for the full alignment
            aln_length_with_gaps = (bp_start+gaps_until_start)
          end
        end
        puts adjusted_domain_frag_pos

        Html.new(frag, html_dir, html_out, aln_aa_html_sub, aln_aa_sub, "#{@dir}/fragments/#{frag}/codeml/", tree_nt, tree_aa, adjusted_domain_frag_pos, project_title, @fasta_file.path, @internal2input_species, @input2internal_species, aln_length_with_gaps, gard_html_file, nucleotide_bias_model, index_html_paths, tex_summary_file.path, tex_summary_file_gapped.path, tex_objects, @refactor, POSEIDON_VERSION, @is_recombination, @timestamp)

        bp_pos += 1
      end
    end


    puts "BUILD MAIN HTML OUTPUT\n"
    ## BUILD MAIN HTML
    html_dir = "#{@dir}/html/"
    html_out = "#{@dir}/html/full_aln/index.html"
    if @root_species
      tree_nt = phylo_nt.tree_rooted
      tree_aa = phylo_aa.tree_rooted
    else
      tree_nt = phylo_nt.tree_corrected
      tree_aa = phylo_aa.tree_corrected
    end

    Html.new('full_aln', html_dir, html_out, aln_aa_html, aln_aa, codeml_full_dir, tree_nt, tree_aa, adjusted_domain_pos, project_title, @fasta_file.path, @internal2input_species, @input2internal_species, 0, gard_html_file, nucleotide_bias_model, index_html_paths, tex_summary_file.path, tex_summary_file_gapped.path, tex_objects, @refactor, POSEIDON_VERSION, @is_recombination, @timestamp)

    puts "BUILD CODEMLHTML\n"
    codeml_html_out = File.open("#{@dir}/html/full_aln/codeml.html",'w')
    CodemlHtml.new('full_aln', codeml_html_out, html_out, tex_objects, frag_names)
    codeml_html_out.close

    puts "BUILD PARAMETER HTML\n"
    parameter_html_out = File.open("#{@dir}/html/full_aln/params.html",'w')
    ParameterHtml.new(parameter_html_out, html_out, frag_names, @timestamp, POSEIDON_VERSION, @parameters)
    parameter_html_out.close

    if fragment_file_paths
      puts "BUILD RECOMBINATION HTML\n"

      html_dir = "#{@dir}/html/"
      html_recomb_out = "#{@dir}/html/full_aln/recomb.html"

      # get full aln tree
      if @root_species
        full_tree_nt = phylo_nt.tree_rooted
        full_tree_aa = phylo_aa.tree_rooted
      else
        full_tree_nt = phylo_nt.tree_corrected
        full_tree_aa = phylo_aa.tree_corrected
      end

      # get fragment trees
      fragment_trees_nt = {}
      fragment_trees_aa = {}
      fragment_file_paths.each do |fragment_file|
        frag = File.dirname(fragment_file).scan(/fragment_[0-9]+/)[0]

        codeml_html_out = File.open("#{@dir}/html/#{frag}/codeml.html",'w')
        CodemlHtml.new(frag, codeml_html_out, html_out, tex_objects, frag_names)
        codeml_html_out.close

        phylo_object = phylo_fragments[fragment_file]
        if @root_species
          fragment_trees_nt[frag] = phylo_object.tree_rooted
          fragment_trees_aa[frag] = phylo_object.tree_rooted.sub('.nt.','.aa.')
        else
          fragment_trees_nt[frag] = phylo_object.tree_corrected
          fragment_trees_aa[frag] = phylo_object.tree_corrected.sub('.nt.','.aa.')
        end
      end

      recomb_html = RecombinationHtml.new(html_dir, html_recomb_out, (gard_html_file+'.adjusted.save'), full_tree_nt, fragment_trees_nt, html_out)
    end


    # archive all results for download link
    Process.fork do
      Dir.chdir(@dir){
        `zip -r #{@timestamp}.zip html/`
      }
    end

    puts "PoSeiDon FINISHED!\n\n"
    @log << "PoSeiDon FINISHED!\n"
    @log.close

  end

  def build_combined_tex(freq2tex, tex_summary_file, insignificant_frags)
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
        frag = File.dirname(tex_file).scan(/fragment_[0-9]+/)[0]


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
    Process.fork do
      Dir.chdir(File.dirname(tex_summary_file.path)) {
        `pdflatex #{tex_summary_file.path}`
        `pdflatex #{tex_summary_file.path}` ## second compilation to adjust the table header width
      }
    end
    Process.waitall
  end

  def build_aln_html_subfile(main_html, frag, bp_start, bp_end, main_aa_aln, gap_start2gap_length)

    puts "#{frag}\t#{bp_start}\t#{bp_end}"
    puts "GAPS: #{gap_start2gap_length}"

    sub_html = File.open("#{@dir}/fragments/#{frag}/aln/#{File.basename(main_html)}",'w')
    sub_aa_aln = File.open("#{@dir}/fragments/#{frag}/aln/#{File.basename(main_aa_aln)}",'w')

    main_nt_aln = main_aa_aln.sub('.aa_ali','.nt_ali')
    sub_nt_aln = File.open("#{@dir}/fragments/#{frag}/aln/#{File.basename(main_aa_aln).sub('.aa_ali','.nt_ali')}",'w')

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

    [sub_html.path, sub_aa_aln.path, aln_length_with_gaps]
  end

  def alignment(fasta_file, aln_out, log)
    cmd = "#{TRANSLATORX_BIN} -i #{fasta_file.path} -p M -o #{aln_out}"
    @parameters['full_aln'] << "#{cmd}\n\n"
    log << "#{cmd}\n"
    `#{cmd}`

    # check the final alignment for N's and if we have some, replace the full codon with ---
    # easiest way: get the X positions from the amino acid files and replace those positions
    aa_aln = Bio::FastaFormat.open("#{aln_out}.aa_ali.fasta")
    nt_aln = Bio::FastaFormat.open("#{aln_out}.nt_ali.fasta")
    html = File.open("#{aln_out}.aa_based_codon_coloured.html",'r')

    @refactor = false
    aa_x_positions = {}
    aa_aln.each do |entry|
      seq = entry.seq
      id = entry.definition
      if seq.include?('X')
        next if seq.count('X') == 1 && seq.reverse.gsub('-','')[0] == 'X'
        @refactor = true
        aa_x_positions[id] = seq.indices('X')
        aa_x_positions[id].delete(seq.length-1) # because we dont want to remove the stop codon also encoded as X
        aa_x_positions.delete(id) if aa_x_positions[id].length == 0
      end
    end
    #aa_x_positions.uniq

    if aa_x_positions.values.join(' ').split(' ').uniq.length > 0
      puts "remove positions #{aa_x_positions} from the amino acid alignment."
      puts "\nPoSeiDon removed positions #{aa_x_positions.values.join(' ').split(' ').uniq} from the amino acid alignment because of ambiguous codons due to N bases in certain species.\n"
	    @mail_notes << "\nPoSeiDon removed positions #{aa_x_positions.values.join(' ').split(' ').uniq} from the amino acid alignment because of ambiguous codons due to N bases in certain species.\n"
    end

    if @refactor
      aa_aln_tmp = File.open("#{aln_out}.aa_ali.fasta.tmp",'w')
      nt_aln_tmp = File.open("#{aln_out}.nt_ali.fasta.tmp",'w')
      html_tmp = File.open("#{aln_out}.aa_based_codon_coloured.html.tmp",'w')

      # refactor aa aln
      Bio::FastaFormat.open("#{aln_out}.aa_ali.fasta").each do |entry|
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

      `cp #{aa_aln.path} #{aa_aln.path}.save`
      `cp #{nt_aln.path} #{nt_aln.path}.save`
      `cp #{html.path} #{html.path}.save`

      `mv #{aa_aln_tmp.path} #{aa_aln.path}`
      `mv #{nt_aln_tmp.path} #{nt_aln.path}`
      `mv #{html_tmp.path} #{html.path}`

    end
  end

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

      fragment_dir = "#{@dir}/fragments/fragment_#{break_counter}/aln/"
      FileUtils.mkdir_p(fragment_dir)

      out = File.open("#{fragment_dir}/#{File.basename(aln)}",'w')
      fragment_file_paths.push(out.path)
      seqs.each do |id, seq|
#        out << ">#{id.chomp}  fragment=#{break_counter} breakpoint=#{(break_position_last.to_i+1).to_s}-#{breakpoint}\n"
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
    break_counter += 1

    fragment_dir = "#{@dir}/fragments/fragment_#{break_counter}/aln/"
    FileUtils.mkdir_p(fragment_dir)

    out = File.open("#{fragment_dir}/#{File.basename(aln)}",'w')
    fragment_file_paths.push(out.path)
    seqs.each do |id, seq|
#      out << ">#{id.chomp}  fragment=#{break_counter} breakpoint=#{(break_position_last.to_i+1).to_s}-#{seq.length}\n"
      out << ">#{id.chomp}\n"
      seq_tmp = seq
      seq_tmp = seq_tmp.sub(seq_splits[id],'')
      out << "#{seq_tmp.scan(/.{1,60}/).join("\n")}\n"
    end
    out.close

    # translate those nt fragments also to aa sequences
    fragment_file_paths.each do |fragment_nt_file|
      fragment_aa_file = File.open(fragment_nt_file.sub('.nt_ali','.aa_ali'),'w')
      in_file = Bio::FastaFormat.open(fragment_nt_file)
      in_file.each do |entry|
        id = entry.definition
        seq = Bio::Sequence::auto(entry.seq)
        aa_seq = seq.translate
        fragment_aa_file << ">#{id.chomp}\n#{aa_seq.chomp}\n"
      end
      fragment_aa_file.close
      in_file.close
    end

    fragment_file_paths
  end

end

class String
  def indices e
    start, result = -1, []
    result << start while start = (self.index e, start + 1)
    result
  end
end



######################################################################################################################################
######################################################################################################################################
## RUN
dir = ARGV[0]
fasta = ARGV[1]

title = ARGV[2].gsub('---',' ')
title = '' if title == 'NA'

root_species = ARGV[3].gsub('---',' ')
if root_species == 'NA'
  root_species = nil
else
  root_species = root_species.split(',')
end

reference_species = ARGV[4].gsub('---',' ')
reference_species = '' if reference_species == 'NA'

kh = ARGV[5]

mail = ARGV[6].gsub('---',' ')
mail = '' if mail == 'NA'

output = ARGV[7]

Main.new(dir, fasta, title, root_species, reference_species, kh, mail, output)
#./main.rb '/mnt/fass2/poseidon-webserver-dev/' '../test_data/bats_mx1.fasta' 'MX1 in bats' '['Rousettus_aegyptiacus', 'Pteropus_alecto', 'Hypsignatus_monstrosus', 'Eidolon_helvum']' '' 'true' 'martin.hoelzer@uni-jena.de' '2019-001'
######################################################################################################################################
######################################################################################################################################


##########################
## TEST CASE
##########################

## normal file, everything nice
fasta = '../test_data/bats_mx1.fasta'
project_title = 'MX1 in bats'
root_species = %w(Rousettus_aegyptiacus Pteropus_alecto Hypsignatus_monstrosus Eidolon_helvum)
query_sequence_name = ''
#Main.new('/mnt/fass2/poseidon-webserver-dev/', fasta, project_title, root_species, query_sequence_name, 'true', 'martin.hoelzer@uni-jena.de', '2019-001')

