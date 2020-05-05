#!/usr/bin/env ruby

# 1) codon colored aln with gaps (orientate on translatex output)
# 2) aa colored aln with gaps
# 3) reorder the sequences in the aln based on the tree (nt and aa), include the trees
# 4) table below alignments with all the models and positive selected positions with p values, color like heatmap by significance
#    positions need to be recalculated to match alignments with the gaps!
# 5) add a position counting (can be implemented as a table with just one row)

require 'bio'
require 'fileutils'

class Html

  def initialize(type, html_dir, out, translatorx_html, aa_aln, codeml_results, nt_tree, aa_tree, domain_pos, title, input_fasta, internal2input_species, input2internal_species, aln_length_with_gaps_adjustor, gard_html_file, nucleotide_bias_model, index_html_paths, tex_summary_file_path, tex_summary_file_path_gapped, tex_objects, refactored_aln, version, is_recomb, logo_png, pipeline, details_summary_master)

    @LOGO_PNG = logo_png
    @PIPELINE = pipeline
    @DETAILS_SUMMARY_MASTER = details_summary_master

    @internal2input_species = internal2input_species
    @input2internal_species = input2internal_species

    FileUtils.mkpath File.dirname(out)
    @out = File.open(out, 'w')

    src_dir = "#{html_dir}/src"
    Dir.mkdir(src_dir) unless Dir.exists?(src_dir)

    # read in the order of species in the nt and aa based tree
    nt_tree_order = read_species_order(nt_tree)
    aa_tree_order = read_species_order(aa_tree)

    # initialize the out file with the nt colored html and initialize parameters for aa alignment
    aa_color_h = {}
    init_html_string = ''
    translatorx_html += '.save' if refactored_aln && !type.include?('frag') #TODO would be nice to replace gaps in the fragment alignment also with the correct N' codons
    codons = init_html(translatorx_html, aa_color_h, init_html_string, nt_tree_order, nt_tree, title, type)
    header_html_string = init_html_string.split('</tr><tr><td></td><td rowspan')[0]
    tmp_string = '<tr><td></td><td rowspan'
    tmp_string << init_html_string.split('</tr><tr><td></td><td rowspan')[1]
    init_html_string = tmp_string

    ## add position row (count codons/aa)
    positions_html_string = '<tr><td></td><td></td>'
    aln_length_with_gaps_adjustor = 1 if type == 'full_aln'
    codons.times do |i|
      positions_html_string << "<td>#{i+aln_length_with_gaps_adjustor}</td>"
    end
    positions_html_string << "</tr>\n"

    # build html aln table for aa aln
    aa_aln_html_string = ''
    aa_aln += '.save' if refactored_aln && !type.include?('frag') #TODO would be nice to replace gaps in the fragment alignment also with the correct N' codons
    gap_positions = aa_html(aa_aln_html_string, aa_aln, aa_color_h, aa_tree_order, aa_tree)
    aa_aln = aa_aln.sub('.save','') if refactored_aln

    ## build a table row with interesting domains
    domain_html_string = '<tr><td></td><td></td>'
    build_domain_html(codons, gap_positions, domain_pos, domain_html_string)

    # build codeml table
    color_gradient = {0.1 => '#DAEAE1', 0.2 => '#C6E9D6', 0.3 => '#C6E9D6', 0.4 => '#B0E5C8', 0.5 => '#B0E5C8', 0.6 => '#9FE5BE', 0.7 => '#9FE5BE', 0.8 => '#4BDD8D', 0.9 => '#28E97F', 0.95 => '#16F179', 0.99 => '#00FF73'}
    codeml_html_string = ''
    result_codeml_files = codeml_html(codeml_html_string, codeml_results, codons, gap_positions, color_gradient, tex_objects, type)

    # copy all files (aln, fasta, tree, ...) to html folder
    html_paths = cp(out, aa_aln, result_codeml_files, nt_tree, input_fasta, gard_html_file, tex_summary_file_path, tex_summary_file_path_gapped, src_dir)

    # copy and refactor the CodeML ctl files, link them in web page
    ctl_refactor(out, aa_aln)

    # add additional content to the left frame
    header_html_string = refactor_framecontent(header_html_string, html_paths, nucleotide_bias_model, type, index_html_paths, version, is_recomb)

    # combine html strings
    @out << header_html_string << domain_html_string << positions_html_string << init_html_string << codeml_html_string << aa_aln_html_string
    @out << '</table></div></div></body></html>'
    @out.close
  end

  def ctl_refactor(out, aa_aln)

    dir_out = "#{File.dirname(out)}/params/codeml_configs"
    FileUtils.mkpath File.dirname(dir_out)

    codeml_dir = File.dirname(aa_aln).sub('/aln','/codeml')

    # copy each ctl file and refactor it
    Dir.glob("#{codeml_dir}/F*/M*").each do |full_path|
      model = full_path.gsub('//','/').split('/').reverse[0]
      freq = full_path.gsub('//','/').split('/').reverse[1]

      out_dir = "#{dir_out}/#{freq}/#{model}"
      FileUtils.mkpath out_dir

      codeml_file = File.open("#{codeml_dir}/#{freq}/#{model}/codeml.ctl",'r')
      codeml_file_refactored = File.open("#{out_dir}/#{freq}_#{model}_codeml.ctl",'w')

      codeml_file.each do |line|
        if line.start_with?('seqfile')
          codeml_file_refactored << "seqfile\t=\t../../../../data/#{File.basename(line.split(' ')[2]).sub('.nogaps.','.nogaps.calc.')}\t\t* sequence data file name\n"
          next
        end
        if line.start_with?('treefile')
          codeml_file_refactored << "treefile\t=\t../../../../data/#{File.basename(line.split(' ')[2])}\t\t* tree structure file name\n"
          next
        end
        if line.start_with?('outfile')
          codeml_file_refactored << "outfile\t=\t#{freq}_#{model}_codeml.mlc\n"
          next
        end
        codeml_file_refactored << line
      end
      codeml_file.close
      codeml_file_refactored.close
    end
  end

  def build_domain_html(codons, gap_positions, domain_pos, domain_html_string)
    col_span = -1
    if domain_pos != 'NA'
      domain_pos.each do |name, range_color_list|
        range_color_list.each do |range_color_a|
          domain_name = name
          domain_color = range_color_a[1]
          domain_start = range_color_a[0].min
          domain_stop = range_color_a[0].max
          col_span = domain_stop - domain_start if range_color_a[2]
          link = "../#{domain_name.downcase}/index.html"
          domain_html_string << "<td style=\"cursor:pointer\" onclick=\"location.href='#{link}'\" bgcolor=\"#{domain_color}\" colspan=\"#{col_span+1}\"><b><font color=\"white\">#{domain_name}</font></b></td>"
        end
      end
    end
    domain_html_string << "</tr>\n\n"
  end

  def refactor_framecontent(header_html_string, html_paths, nucleotide_bias_model, type, index_html_paths, version, is_recomb)
    data_html = "</br>\n<details><summary><b>Input</b></summary><ul><li><a href=\"#{html_paths[:input_fasta]}\">Fasta (nt)</a></li><li><a href=\"#{html_paths[:input_fasta_aa]}\">Fasta (aa)</a></li></ul></details>"

    data_html << "\n<details><summary><b>Alignment</b></summary><ul><li><b>nt</b>:<br><a href=\"#{html_paths[:aln_nt]}\">full</a>, <a href=\"#{html_paths[:aln_nt_nogap]}\">no gaps</a></li><li><b>aa</b>:<br><a href=\"#{html_paths[:aln_aa]}\">full</a>, <a href=\"#{html_paths[:aln_aa_nogap]}\">no gaps</a></li></ul></details>"

    data_html << "\n<details><summary><b>Tree</b></summary><ul><li><b>nt</b>:<br> <a href=\"#{html_paths[:tree_nt_newick]}\">newick</a>, <a href=\"#{html_paths[:tree_nt_scaled_pdf]}\">pdf</a>, <a href=\"#{html_paths[:tree_nt_scaled_svg]}\">svg</a></li><li><b>aa</b>:<br> <a href=\"#{html_paths[:tree_aa_newick]}\">newick</a>, <a href=\"#{html_paths[:tree_aa_scaled_pdf]}\">pdf</a>, <a href=\"#{html_paths[:tree_aa_scaled_svg]}\">svg</a></li></ul></details>"

    data_html << "\n<details><summary><b>Codeml</b></summary><ul><li><a href=\"codeml.html\">Go to overview</a></li></ul></details>"
#    data_html << "\n<details><summary><b>Codeml</b></summary><ul>"
#    html_paths[:codeml].each do |codeml|
#      ns_sites = File.dirname(codeml).split('/').reverse[0]
#      freq = File.dirname(codeml).split('/').reverse[1]
#      data_html << "<li>#{freq}, #{ns_sites}"
#      data_html << "<ul><li><a href=\"#{codeml}\">results</a> | <a href=\"params/codeml_configs/#{freq}/#{ns_sites}/#{freq}_#{ns_sites}_codeml.ctl\">config</a></li></ul></li>"
#    end
#    data_html << '</ul>'
#    data_html << '</ul></details>'
    
    if type.include?('fragment_')
      data_html << "\n<details open=\"open\"><summary><b>Recombination</b></summary><ul><li>Model: <a href=\"#{html_paths[:model]}\">#{nucleotide_bias_model}</a></li>"
      data_html << "<li><a href=\"../full_aln/recomb.html\">Tree topologies</a></li>" if is_recomb
      data_html << "<li><a href=\"#{html_paths[:gard]}\">GARD</a></li>"
    else
      data_html << "\n<details><summary><b>Recombination</b></summary><ul><li>Model: <a href=\"#{html_paths[:model]}\">#{nucleotide_bias_model}</a></li>"
      data_html << "<li><a href=\"recomb.html\">Tree topologies</a></li>" if is_recomb
      data_html << "<li><a href=\"#{html_paths[:gard]}\">GARD</a></li>"
    end
    # loop over fragments and link them here
    if index_html_paths.keys.size > 1 # then we have fragments and not only the full aln html
      data_html << '<ul>'
      index_html_paths.each do |frag_type, html_path|
        unless frag_type == type
          data_html << "<li><a href=\"#{html_path}\">#{frag_type}</a></li>"
        end
      end
      data_html << '</ul>'
    end
    data_html << '</ul></details>'

    data_html << "\n<details open=\"open\"><summary><b>Significant sites</b></summary><ul><li><b>no gaps</b>:<br><a href=\"#{html_paths[:tex_tex]}\">tex</a>, <a href=\"#{html_paths[:tex_pdf]}\">pdf</a></li>"
    data_html << "\n<li><b>with gaps</b>:<br><a href=\"#{html_paths[:tex_tex_gapped]}\">tex</a>, <a href=\"#{html_paths[:tex_pdf_gapped]}\">pdf</a></li></ul></details>"
    if type.include?('fragment_')
      data_html << "\n<a href=\"../full_aln/index.html\">go back to full alignment</a>"
    end

    if type.include?('fragment_')
      data_html << "\n<div class=\"bottom\"><p><a href=\"../full_aln/params.html\"<b>Parameters</b></a></p>\n"
    else
      data_html << "\n<div class=\"bottom\"><p><a href=\"params.html\"<b>Parameters</b></a></p>\n"
    end
    data_html << "<p>PoSeiDon v#{version}</p></div>"

    # insert the logo
    #data_html << "\n<a target=\"_blank\" href=\"http://www.rna.uni-jena.de/en/poseidon\"><img src=\"#{File.basename(@LOGO_PNG)}\" width=\"120px\"/></a>"

    header1 = header_html_string.split('</table>')[0].gsub('100px; /','145px; /')
    header1 = header1.sub('145px; /','170px; /')
    header2 = header_html_string.split('</table>')[1]
    header1 + '</table>' + data_html + "\n\n" + header2
  end

  def cp(out, aa_aln, codeml_results, nt_tree, input_fasta, gard_html_file, tex_summary_file_path, tex_summary_file_path_gapped, src_dir)

    file_out = "#{File.dirname(out)}/data"
    Dir.mkdir(file_out) unless Dir.exists?(file_out)

    html_paths = {}

    # the logo & pipeline
    `cp #{@LOGO_PNG} #{src_dir}/poseidon_logo.png` unless File.exists?("#{src_dir}/poseidon_logo.png")
    `cp #{@PIPELINE}.* #{src_dir}/` unless File.exists?("#{src_dir}/#{@PIPELINE}.png")

    # the details-summary-master folder
    `cp -r #{@DETAILS_SUMMARY_MASTER} #{src_dir}/` unless Dir.exists?("#{src_dir}/#{File.basename(@DETAILS_SUMMARY_MASTER)}")

    # tex summary files (with gaps and without)
    `cp #{tex_summary_file_path} #{file_out}/`; html_paths[:tex_tex] = "data/#{File.basename(tex_summary_file_path)}"
    `cp #{tex_summary_file_path.sub('.tex','.pdf')} #{file_out}/`; html_paths[:tex_pdf] = "data/#{File.basename(tex_summary_file_path).sub('.tex','.pdf')}"
    `cp #{tex_summary_file_path_gapped} #{file_out}/`; html_paths[:tex_tex_gapped] = "data/#{File.basename(tex_summary_file_path_gapped)}"
    `cp #{tex_summary_file_path_gapped.sub('.tex','.pdf')} #{file_out}/`; html_paths[:tex_pdf_gapped] = "data/#{File.basename(tex_summary_file_path_gapped).sub('.tex','.pdf')}"

    # get all input fastas
    `cp #{input_fasta} #{file_out}`; refactor_species_names("#{file_out}/#{File.basename(input_fasta)}", @internal2input_species); html_paths[:input_fasta] = "data/#{File.basename(input_fasta)}";
    if File.exists?(aa_aln.sub('_ali','seqs'))
      `cp #{aa_aln.sub('_ali','seqs')} #{file_out}`; refactor_species_names("#{file_out}/#{File.basename(aa_aln.sub('_ali','seqs'))}", @internal2input_species); html_paths[:input_fasta_aa] = "data/#{File.basename(aa_aln.sub('_ali','seqs'))}"
    end

    # get all aln
    `cp #{aa_aln.sub('.aa','.nt')} #{file_out}`; refactor_species_names("#{file_out}/#{File.basename(aa_aln.sub('.aa','.nt'))}", @internal2input_species); html_paths[:aln_nt] = "data/#{File.basename(aa_aln.sub('.aa','.nt'))}"
    `cp #{aa_aln} #{file_out}`; refactor_species_names("#{file_out}/#{File.basename(aa_aln)}", @internal2input_species); html_paths[:aln_aa] = "data/#{File.basename(aa_aln)}"
    `cp #{aa_aln.sub('.aa_ali.fasta','.nt_ali.nogaps.fasta')} #{file_out}`; refactor_species_names("#{file_out}/#{File.basename(aa_aln.sub('.aa_ali.fasta','.nt_ali.nogaps.fasta'))}", @internal2input_species); html_paths[:aln_nt_nogap] = "data/#{File.basename(aa_aln.sub('.aa_ali.fasta','.nt_ali.nogaps.fasta'))}"
    `cp #{aa_aln.sub('.aa_ali.fasta','.nt_ali.nogaps.fasta')} #{file_out}/#{File.basename(aa_aln).sub('.aa_ali.fasta','.nt_ali.nogaps.calc.fasta')}` # copy again for reproducable calculations
    `cp #{aa_aln.sub('_ali.fasta','_ali.nogaps.fasta')} #{file_out}`; refactor_species_names("#{file_out}/#{File.basename(aa_aln.sub('_ali.fasta','_ali.nogaps.fasta'))}", @internal2input_species); html_paths[:aln_aa_nogap] = "data/#{File.basename(aa_aln.sub('_ali.fasta','_ali.nogaps.fasta'))}"

    #TODO  get all trees and build PDF, SVG and PNG files on original species names

    if File.basename(nt_tree).include?('.rooted')
      `cp #{nt_tree.sub('.corrected.rooted','')} #{file_out}/#{File.basename(nt_tree).sub('.corrected.rooted','.newick')}`; html_paths[:tree_nt_newick] = "data/#{File.basename(nt_tree.sub('.corrected.rooted','.newick'))}"
      `cp #{nt_tree.sub('.corrected.rooted','.rooted.pdf')} #{file_out}`; html_paths[:tree_nt_pdf] = "data/#{File.basename(nt_tree.sub('.corrected.rooted','.pdf'))}"
      `cp #{nt_tree.sub('.corrected.rooted','.rooted.svg')} #{file_out}`; html_paths[:tree_nt_svg] = "data/#{File.basename(nt_tree.sub('.corrected.rooted','.svg'))}"
      `cp #{nt_tree.sub('.corrected.rooted','.rooted.png')} #{file_out}`
      `cp #{nt_tree.sub('.corrected.rooted','.unshod.tree')} #{file_out}`
      `cp #{nt_tree.sub('.nt','.aa').sub('.corrected.rooted','')} #{file_out}/#{File.basename(nt_tree).sub('.nt','.aa').sub('.corrected.rooted','.newick')}`; html_paths[:tree_aa_newick] = "data/#{File.basename(nt_tree.sub('.nt','.aa').sub('.corrected.rooted','.newick'))}"
      `cp #{nt_tree.sub('.nt','.aa').sub('.corrected.rooted','.rooted.pdf')} #{file_out}`; html_paths[:tree_aa_pdf] = "data/#{File.basename(nt_tree.sub('.nt','.aa').sub('.corrected.rooted','.pdf'))}"
      `cp #{nt_tree.sub('.nt','.aa').sub('.corrected.rooted','.rooted.svg')} #{file_out}`; html_paths[:tree_aa_svg] = "data/#{File.basename(nt_tree.sub('.nt','.aa').sub('.corrected.rooted','.svg'))}"
      `cp #{nt_tree.sub('.nt','.aa').sub('.corrected.rooted','.rooted.png')} #{file_out}`
      # get also the trees with the scale to link them as downloadable svgs and pdfs
      `cp #{nt_tree.sub('.corrected.rooted','.rooted.scale.pdf')} #{file_out}`; html_paths[:tree_nt_scaled_pdf] = "data/#{File.basename(nt_tree.sub('.corrected.rooted','.rooted.scale.pdf'))}"
      `cp #{nt_tree.sub('.corrected.rooted','.rooted.scale.svg')} #{file_out}`; html_paths[:tree_nt_scaled_svg] = "data/#{File.basename(nt_tree.sub('.corrected.rooted','.rooted.scale.svg'))}"
      `cp #{nt_tree.sub('.nt','.aa').sub('.corrected.rooted','.rooted.scale.pdf')} #{file_out}`; html_paths[:tree_aa_scaled_pdf] = "data/#{File.basename(nt_tree.sub('.nt','.aa').sub('.corrected.rooted','.rooted.scale.pdf'))}"
      `cp #{nt_tree.sub('.nt','.aa').sub('.corrected.rooted','.rooted.scale.svg')} #{file_out}`; html_paths[:tree_aa_scaled_svg] = "data/#{File.basename(nt_tree.sub('.nt','.aa').sub('.corrected.rooted','.rooted.scale.svg'))}"
    else
      `cp #{nt_tree.sub('.corrected','')} #{file_out}/#{File.basename(nt_tree).sub('.corrected','.newick')}`; html_paths[:tree_nt_newick] = "data/#{File.basename(nt_tree.sub('.corrected','.newick'))}"
      `cp #{nt_tree.sub('.corrected','.pdf')} #{file_out}`; html_paths[:tree_nt_pdf] = "data/#{File.basename(nt_tree.sub('.corrected','.pdf'))}"
      `cp #{nt_tree.sub('.corrected','.svg')} #{file_out}`; html_paths[:tree_nt_svg] = "data/#{File.basename(nt_tree.sub('.corrected','.svg'))}"
      `cp #{nt_tree.sub('.corrected','.png')} #{file_out}`
      `cp #{nt_tree.sub('.corrected.rooted','.unshod.tree')} #{file_out}`
      `cp #{nt_tree.sub('.nt','.aa').sub('.corrected','')} #{file_out}/#{File.basename(nt_tree).sub('.nt','.aa').sub('.corrected','.newick')}`; html_paths[:tree_aa_newick] = "data/#{File.basename(nt_tree.sub('.nt','.aa').sub('.corrected','.newick'))}"
      `cp #{nt_tree.sub('.nt','.aa').sub('.corrected','.pdf')} #{file_out}`; html_paths[:tree_aa_pdf] = "data/#{File.basename(nt_tree.sub('.nt','.aa').sub('.corrected','.pdf'))}"
      `cp #{nt_tree.sub('.nt','.aa').sub('.corrected','.svg')} #{file_out}`; html_paths[:tree_aa_svg] = "data/#{File.basename(nt_tree.sub('.nt','.aa').sub('.corrected','.svg'))}"
      `cp #{nt_tree.sub('.nt','.aa').sub('.corrected','.png')} #{file_out}`
      # get also the trees with the scale to link them as downloadable svgs and pdfs
      `cp #{nt_tree.sub('.corrected','.scale.pdf')} #{file_out}`; html_paths[:tree_nt_scaled_pdf] = "data/#{File.basename(nt_tree.sub('.corrected','.scale.pdf'))}"
      `cp #{nt_tree.sub('.corrected','.scale.svg')} #{file_out}`; html_paths[:tree_nt_scaled_svg] = "data/#{File.basename(nt_tree.sub('.corrected','.scale.svg'))}"
      `cp #{nt_tree.sub('.nt','.aa').sub('.corrected','.scale.pdf')} #{file_out}`; html_paths[:tree_aa_scaled_pdf] = "data/#{File.basename(nt_tree.sub('.nt','.aa').sub('.corrected','.scale.pdf'))}"
      `cp #{nt_tree.sub('.nt','.aa').sub('.corrected','.scale.svg')} #{file_out}`; html_paths[:tree_aa_scaled_svg] = "data/#{File.basename(nt_tree.sub('.nt','.aa').sub('.corrected','.scale.svg'))}"
    end

    # get all codeml result files
    html_paths[:codeml] = []
    codeml_results.each do |codeml|
      ns_sites = File.dirname(codeml).split('/').reverse[0]
      freq = File.dirname(codeml).split('/').reverse[1]

      `mkdir -p #{file_out}/codeml/#{freq}/#{ns_sites}; cp #{codeml} #{file_out}/codeml/#{freq}/#{ns_sites}/`
      html_paths[:codeml].push("data/codeml/#{freq}/#{ns_sites}/#{File.basename(codeml)}")
    end

    # get GARD html file and model selection log file
    html_paths[:gard] = "data/#{File.basename(gard_html_file)}"
    `cp #{gard_html_file} #{file_out}/#{File.basename(gard_html_file)}`
    html_paths[:model] = "data/model.log"
    `cp #{File.dirname(gard_html_file).sub('/gard','/model_selection')}/model.log #{file_out}/model.log`

    html_paths
  end

  def read_species_order(tree)
    input = Bio::FlatFile.open(Bio::Newick, tree).next_entry.tree
    input.leaves
  end

  def codeml_html(out, codeml_results, codons_counter, gap_positions, color_gradient, tex_objects, type)

    freqs = ['F3X4', 'F1X4', 'F61']
    models = ['M0', 'M1a', 'M2a', 'M7', 'M8', 'M8a']

    result_files = []
    # for each leave folder, check if in the outfile are positive selected calculations are included and if so: build one table row
    freqs.each do |freq|
      models.each do |model|
        ns_sites_folder = "#{codeml_results}/#{freq}/#{model}"
        ns_sites = model
        codon_freq = freq
        # check if outfile includes positive selection calcs
        pos_selec_calcs = false
        pos_h = {}
        File.open("#{ns_sites_folder}/codeml.mlc",'r').each do |line|
          if line.include?('Bayes Empirical Bayes (BEB) analysis')
            pos_selec_calcs = true
          end
          result_files.push("#{ns_sites_folder}/codeml.mlc") unless result_files.include?("#{ns_sites_folder}/codeml.mlc")
          if pos_selec_calcs
#            result_files.push("#{ns_sites_folder}/codeml.mlc") unless result_files.include?("#{ns_sites_folder}/}.mlc")
            split = line.split(' ')
            if split[0].to_i > 0
              pos = split[0].to_i
              aa = split[1]
              probability = split[2]

              # adjust position because we have gaps in the alignment
              old_adjuster = 0
              adjuster = 0
              gap_positions.each do |gap_pos|
                adjuster += 1 if gap_pos <= pos
              end
              #pos = pos + adjuster
              while old_adjuster != adjuster
                old_adjuster = adjuster
                adjuster = adjust_pos(pos, old_adjuster, gap_positions)
              end
              pos = pos + adjuster - 1

              pos_h[pos] = [aa, probability]
            end
          end
        end

        if pos_h.keys.size > 0 || ns_sites == 'M2a' || ns_sites == 'M8'

        case ns_sites
          when 'M2a'
            null_model_type = 'M1a'
          when 'M8'
            null_model_type = 'M7/M8a'
          else
            null_model_type = ''
        end

        lrt = ''; pvalue = ''; lnL0 = ''; lnL1 = ''; omega_percent = ''; omega_average = ''
        m0 = ''; m1 = ''
        lrt_2 = ''; pvalue_2 = ''; lnL0_2 = ''; m0_2 = ''
        tex_objects[type].each do |tex_object|
          next unless tex_object.include?('.gaps.')
          # here we have now the gapped und ungapped tex files!
          # tex_object="tex/codeml_F3X4.all.M1a_vs_M2a.tex"
            #new_freq = File.dirname(tex_object.output_file.path).reverse.split('/')[0].reverse
            new_freq = File.basename(tex_object, '.tex').sub('codeml_','').split('.')[0]

            if new_freq == codon_freq
              if ns_sites == 'M2a'
                f = File.open("#{codon_freq}.lrt",'r')
                f.each do |l|
                  s = l.split("\t")
                  lrt = s[1].chomp if s[0] == 'm2a_lrt'
                  pvalue = s[1].gsub('$','').sub('<','').to_f if s[0] == 'm2a_pvalue'
                  lnL0 = s[1].chomp if s[0] == 'm2a_lnL0'
                  lnL1 = s[1].chomp if s[0] == 'm2a_lnL1'
                  omega_percent = s[1].chomp if s[0] == 'm2a_omega_percent'
                  omega_average = s[1].chomp if s[0] == 'm2a_omega_average'
                  m0 = 'M1a'
                  m1 = 'M2a'
                end
                f.close
                #lrt = tex_object.significance_m2a.lrt
                #pvalue = tex_object.significance_m2a.pvalue.gsub('$','').sub('<','').to_f
                #lnL0 = tex_object.significance_m2a.lnL0
                #lnL1 = tex_object.significance_m2a.lnL1
                #omega_percent = tex_object.significance_m2a.omega_percent
                #omega_average = tex_object.significance_m2a.average_omega
                #m0 = 'M1a'
                #m1 = 'M2a'
              end
              if ns_sites == 'M8'
                f = File.open("#{codon_freq}.lrt",'r')
                f.each do |l|
                  s = l.split("\t")
                  lrt = s[1].chomp if s[0] == 'm8_lrt'
                  lrt_2 = s[1].chomp if s[0] == 'm8_lrt_2'
                  pvalue = s[1].gsub('$','').sub('<','').to_f if s[0] == 'm8_pvalue'
                  pvalue_2 = s[1].gsub('$','').sub('<','').to_f if s[0] == 'm8_pvalue_2'
                  lnL0 = s[1].chomp if s[0] == 'm8_lnL0'
                  lnL0_2 = s[1].chomp if s[0] == 'm8_lnL0_2'
                  lnL1 = s[1].chomp if s[0] == 'm8_lnL1'
                  omega_percent = s[1].chomp if s[0] == 'm8_omega_percent'
                  omega_average = s[1].chomp if s[0] == 'm8_omega_average'
                  m0 = 'M7'
                  m0_2 = 'M8a'
                  m1 = 'M8'
                end
                f.close
                #lrt = tex_object.significance_m8.lrt
                #lrt_2 = tex_object.significance_m8a.lrt
                #pvalue = tex_object.significance_m8.pvalue.gsub('$','').sub('<','').to_f
                #pvalue_2 = tex_object.significance_m8a.pvalue.gsub('$','').sub('<','').to_f
                #lnL0 = tex_object.significance_m8.lnL0
                #lnL0_2 = tex_object.significance_m8a.lnL0
                #lnL1 = tex_object.significance_m8.lnL1
                #omega_percent = tex_object.significance_m8.omega_percent
                #omega_average = tex_object.significance_m8.average_omega
                #m0 = 'M7'
                #m0_2 = 'M8a'
                #m1 = 'M8'
            end
          end
        end
        puts "#{ns_sites}\t#{codon_freq}\t#{lrt}\t#{pvalue}\t#{lnL0}\t#{lnL1}\t#{omega_percent}\t#{omega_average}"

        color = '#A9F5A9'
        comparison = '='
        comparison_2 = '='
        text = 'LRT is significant'
        if pvalue < 0.1
          comparison = '<' if pvalue == 0.001
        else
          color = '#FA5858'
          text = 'LRT not significant'
          omega_percent = 'NA'
          omega_average = 'NA'
        end

        # check again because of second null model M8a
        if ns_sites == 'M8'
          if pvalue_2 < 0.1
            comparison_2 = '<' if pvalue_2 == 0.001
            if text == 'LRT not significant'
              text = 'LRT is mixed!'
              color = '#B45F04'
            end
          else
            if text == 'LRT is significant'
              text = 'LRT is mixed!'
              color = '#B45F04'
            end
          end
        end

        if ns_sites == 'M2a'
          out << "<tr><td></td><td class=\"alnleft\" bgcolor=\"#{color}\">
<details><summary>#{codon_freq}, #{null_model_type} vs #{ns_sites} | #{text}</summary>
<br>lnL_<sub>#{m0}</sub> = #{lnL0} <br> lnL_<sub>#{m1}</sub> = #{lnL1} <br> LRT: &chi;<sup>2</sup> = #{lrt} <br> % sites with &omega; > 1 = #{omega_percent} <br> avg(&omega;) = #{omega_average} <br> p #{comparison} #{pvalue}</details></td>\n"
        else
          # show combined results of M7 and M8a vs M8 test
          out << "<tr><td></td><td class=\"alnleft\" bgcolor=\"#{color}\">
<details><summary>#{codon_freq}, #{null_model_type} vs #{ns_sites} | #{text}</summary>
<br>lnL_<sub>#{m0}</sub> = #{lnL0} / lnL_<sub>#{m0_2}</sub> = #{lnL0_2} <br> lnL_<sub>#{m1}</sub> = #{lnL1} <br> LRT: &chi;<sup>2</sup> = #{lrt} / LRT: &chi;<sup>2</sup> = #{lrt_2} <br> % sites with &omega; > 1 = #{omega_percent} <br> avg(&omega;) = #{omega_average} <br> p #{comparison} #{pvalue} / p #{comparison_2} #{pvalue_2}</details></td>\n"
        end

          # we found values in this calculation, fill table row
          codons_counter.times do |i|
            if pos_h[i]
              probability = pos_h[i][1].gsub('*','')
              html_color = '#FFFFFF'
              color_gradient.each do |cutoff, color|
                html_color = color if probability.to_f >= cutoff
              end
              if ns_sites == 'M2a'
                out << "<td bgcolor=\"#{html_color}\"><span title=\"#{codon_freq}, #{null_model_type} vs #{ns_sites} | #{text} (p #{comparison} #{pvalue})\">#{probability}</span></td>\n"
              else
                out << "<td bgcolor=\"#{html_color}\"><span title=\"#{codon_freq}, #{null_model_type} vs #{ns_sites} | #{text} (p #{comparison} #{pvalue} / p #{comparison_2} #{pvalue_2})\">#{probability}</span></td>\n"
              end
            else
              out << "<td></td>\n"
            end
          end
          out << "</tr>\n"
        end
      end
    end
    result_files
  end

  def adjust_pos(pos, old_adjuster, gap_positions)
    adjuster = 0
    gap_positions.each do |gap_pos|
      adjuster += 1 if gap_pos <= (pos+old_adjuster)
    end
    adjuster
  end

  def refactor_species_names(file, hash)
    old = File.open(file,'r')
    tmp = File.open("#{file}.tmp",'w')
    old_string = ''
    old.each do |l|
      old_string += l
    end
    hash.each do |species_name1, species_name2|
      old_string = old_string.sub(species_name1,species_name2)
    end
    tmp << old_string
    tmp.close
    old.close
    `mv #{tmp.path} #{old.path}`
  end

  def aa_html(out, aa_aln, aa_color_h, aa_tree_order, aa_tree)
    gap_positions = []
    seqs = {}
    Bio::FastaFormat.open(aa_aln).each do |entry|
      split = entry.definition.split('_')
      search_species = split.join(' ')
      seqs[search_species] = entry.seq
    end

    # reorder species by aa tree
    seqs_ordered = {}
    aa_tree_order.each do |species|
#      puts species.name
#      puts seqs[species.name]
      seqs_ordered[species.name] = seqs[species.name]
    end

    picture_is_placed = false

    seqs_ordered.each do |id, seq|
      if picture_is_placed
        out << "</tr><tr><td></td>\n"
      else
        if aa_tree.include?('.rooted')
          out << "<tr><td></td><td rowspan=\"#{aa_tree_order.size}\"><a target=\"_blank\" href=\"data/#{File.basename(aa_tree).sub('.corrected.rooted','.rooted.scale.svg')}\"><img src=\"data/#{File.basename(aa_tree).sub('.corrected.rooted','.rooted.png')}\" /></a></td>"
        else
          out << "<tr><td></td><td rowspan=\"#{aa_tree_order.size}\"><a target=\"_blank\" href=\"data/#{File.basename(aa_tree).sub('.corrected','.scale.svg')}\"><img src=\"data/#{File.basename(aa_tree).sub('.corrected','.png')}\" /></a></td>"
        end
        picture_is_placed = true
      end

      c = 0
      seq.scan(/.{0,1}/).each do |aa|
        species_name = @internal2input_species[id.gsub(' ','_').chomp]
        out << "<td style=\"height:23px\" bgcolor=\"#{aa_color_h[aa.intern]}\"><span title=\"#{species_name}\">#{aa}</span></td>\n"
        c += 1
        gap_positions.push(c) if aa == '-'
      end
    end
    out << "</tr>\n\n"
    gap_positions.uniq.sort
  end

  def init_html(translatorx_html, aa_color_h, init_html_string, tree_order, tree_path, title, type)
    read_header = true
    read_aa_colors = false
    codons = 0
    taxa = 0
    count_codons = false
    picture_is_placed = false

    species_row = {}; read_new_species_row = false
    header = '' # all before the species rows
    footer = '' # all after the species rows

    species = nil
    File.open(translatorx_html).each do |html_line|

      count_codons = true if html_line.include?('id="maincontent"')
      codons += 1 if html_line.include?('bgcolor') && count_codons
      #count_codons = false if html_line.include?('TAA') || html_line.include?('TAG') || html_line.include?('TGA')
      taxa += 1 if html_line.include?('</tr><tr><td>') && count_codons

      if html_line.include?('</style>')
        header << "\t.innertube img {\n\t\tposition: relative;\n\t\tleft: -10;\n\t\ttop: -5;\n\t}"
        header << "      #alnTable td
    {
    	text-align:center;
    	vertical-align:middle;
		}\n
    /* unvisited link */
    a:link {
      color: white;
    }

    /* visited link */
    a:visited {
      color: white;
    }

    details {
      cursor: pointer;
    }

    a.nounderline {
			text-decoration:none;
    }

    .alnleft { text-align: left!important; }

    div.bottom {
        position: relative;
        top: 100px;
    }

    ul {
      margin:12;
      padding: 0;
    }"
      end

      if html_line.start_with?('</tr><tr><td>')
        read_new_species_row = true
        species = html_line.split('<td>')[1].split('</td>')[0].gsub('_',' ')
        species_row[species] = "</tr><tr><td></td>\n"
      else
        if html_line.include?('</table></div></div></body></html>')
          footer << html_line.sub('</table></div></div></body></html>','').sub('table border=0><tr>','table border=0 id="alnTable"><tr></tr>')
        else
          if html_line.include?('bgcolor') && read_new_species_row
            species_row[species] << html_line.sub('">',"\" style=\"height:23px\"><span title=\"#{@internal2input_species[species.gsub(' ','_').chomp]}\">").sub('</td>','</span></td>')
          else
            if html_line.include?('<div id="framecontent"><div class="innertube">')
              subtitle = 'full alignment'
              if type.include?('fragment')
                subtitle = type
              end
              header << "<div id=\"framecontent\"><div class=\"innertube\">\n<h3>#{title}</h3><h4>#{subtitle}</h4>\n"
            else
              header << html_line.sub('table border=0><tr>','table border=0 id="alnTable"><tr></tr>').sub('AA-coloured codon alignment',"#{title}")
            end
          end
        end
      end

      if read_aa_colors
        unless html_line.start_with?('</tr>')
          aa = html_line.split('>')[1].sub('</td','').intern
          color = html_line.scan(/#\w+/)[0]
          aa_color_h[aa] = color
        end
        read_aa_colors = false if html_line.include?('</table>')
      end

      if html_line.include?('<table><tr>')
        read_aa_colors = true
      end
    end

    # reorder species aln rows
    ordered_aln_string = ''
    picture_is_placed = true
    tree_order.each do |species|
      if picture_is_placed
        picture_is_placed = false
        #puts species_row.keys
        if tree_path.include?('.rooted')
          ordered_aln_string << species_row[species.name].sub('</tr><tr><td></td>',"</tr><tr><td></td><td rowspan=\"#{tree_order.size}\"><a target=\"_blank\" href=\"data/#{File.basename(tree_path).sub('.corrected.rooted','.rooted.scale.svg')}\"><img src=\"data/#{File.basename(tree_path).sub('.corrected.rooted','.rooted.png')}\" /></a></td>")
        else
          ordered_aln_string << species_row[species.name].sub('</tr><tr><td></td>',"</tr><tr><td></td><td rowspan=\"#{tree_order.size}\"><a target=\"_blank\" href=\"data/#{File.basename(tree_path).sub('.corrected','.scale.svg')}\"><img src=\"data/#{File.basename(tree_path).sub('.corrected','.png')}\" /></a></td>")
        end
      else
        ordered_aln_string << species_row[species.name]
      end
    end

    header = header.sub('overflow: hidden; /*','overflow: scroll; /*')
    #add logo
    header = header.split('<h3>')[0] << "\n<a target=\"_blank\" href=\"http://www.rna.uni-jena.de/en/poseidon\"><img src=\"../src/#{File.basename(@LOGO_PNG)}\" width=\"120px\"/></a>\n" << '<h3>' << header.split('<h3>')[1]
    #add details-summary-master sources
    header1 = header.split('</style>')[0]
    header2 = header.split('</style>')[1]
    header = header1 << '</style>' << "\n\n<script type=\"text/javascript\" src=\"../src/details-shim-master/details-shim.min.js\">\n</script><link rel=\"stylesheet\" type=\"text/css\" href=\"../src/details-shim-master/details-shim.min.css\">\n\n" << header2

    init_html_string << header << ordered_aln_string << footer
    codons / taxa
  end

end

##################################################################
## INPUT

type = ARGV[0]
html_dir = ARGV[1]
out = ARGV[2]
translatorx_html = ARGV[3]
aa_aln = ARGV[4]
codeml_results = ARGV[5] # only include the ctl atm
nt_tree = ARGV[6]
aa_tree = ARGV[7]
domain_pos = ARGV[8] # should be a hash, 'NA' if not set 
title = ARGV[9]
input_fasta = ARGV[10]

internal2input_species_tsv = File.open(ARGV[11],'r')
input2internal_species_tsv = File.open(ARGV[12],'r')
internal2input_species = {}
internal2input_species_tsv.each do |line|
  internal2input_species[line.split("\t")[0]] = line.split("\t")[1].chomp 
end
internal2input_species_tsv.close
input2internal_species = {}
input2internal_species_tsv.each do |line|
  input2internal_species[line.split("\t")[0]] = line.split("\t")[1].chomp 
end
input2internal_species_tsv.close

aln_length_with_gaps_adjustor = ARGV[13].to_i

gard_html_file = ARGV[14]
nucleotide_bias_model = ARGV[15]

# hash holding the paths to the full and fragment index.html files
# generate this based on the html files that we have, for now only for full aln
index_html_paths = {} 
index_html_paths['full_aln'] = '../full_aln/index.html'

tex_summary_file_path = ARGV[16] 
tex_summary_file_path_gapped = ARGV[17]

# TODO: Fragments!!
tex_dir = ARGV[18]
tex_objects = {'full_aln' => []}
Dir.glob("#{tex_dir}/*.tex").each do |tex_full_aln|
  tex_objects['full_aln'].push(tex_full_aln)
end

refactored_aln = ARGV[19].to_s.downcase == "true"

version = ARGV[20]

is_recomb = ARGV[21].to_s.downcase == "true"

workflow_projectdir = ARGV[22]
logo = "#{workflow_projectdir}/images/poseidon_logo.png"
pipeline = "#{workflow_projectdir}/images/pipeline_landscape"
details = "#{workflow_projectdir}/src/details-shim-master"

Html.new(type, html_dir, out, translatorx_html, aa_aln, codeml_results, nt_tree, aa_tree, domain_pos, title, input_fasta, internal2input_species, input2internal_species, aln_length_with_gaps_adjustor, gard_html_file, nucleotide_bias_model, index_html_paths, tex_summary_file_path, tex_summary_file_path_gapped, tex_objects, refactored_aln, version, is_recomb, logo, pipeline, details)

######################
## TEST: MX1 bat species, Jonas cloned ones and additional bats, only one allel per species (see manuscript table)
######################
#title = 'MX1 in bat species'
#html_out = '/home/hoelzer/projects/mx_bat_georg/FINAL_CALCS_120516/FINAL_DATA/html/index.html'
#translatorx_html_output = '/home/hoelzer/projects/mx_bat_georg/FINAL_CALCS_120516/aln/bats_all/bats_all_sorted.aa_based_codon_coloured.html'
#translatorx_aa_aln = '/home/hoelzer/projects/mx_bat_georg/FINAL_CALCS_120516/aln/bats_all/bats_all_sorted.aa_ali.fasta'
#codeml_results = '/home/hoelzer/projects/mx_bat_georg/FINAL_CALCS_120516/codeml_for_html/bats_all_sorted/'
#nt_tree = '/home/hoelzer/projects/mx_bat_georg/FINAL_CALCS_120516/phylo/bats_all_nogaps/RAxML_bipartitionsBranchLabels.nt.corrected'
#aa_tree = '/home/hoelzer/projects/mx_bat_georg/FINAL_CALCS_120516/phylo/bats_all_nogaps/RAxML_bipartitionsBranchLabels.nt.corrected'
#input_fasta = '/home/hoelzer/projects/mx_bat_georg/FINAL_CALCS_120516/data/bats_all.fasta.sorted'
#domain_pos = {:BSE => [[(18..39),:red,true],[(312..337),:red,true],[(571..600),:red,true]], 'G-DOMAIN' => [[(40..311), :orange,true]], :STALK => [[(338..501), :blue,true], [(523..570), :blue,true]], 'LOOP L4' => [[(489..532), :green,true]]}
#Html.new(html_out, translatorx_html_output, translatorx_aa_aln, codeml_results, nt_tree, aa_tree, domain_pos, title, input_fasta)
