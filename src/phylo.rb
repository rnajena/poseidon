#!/home/hoelzer/local/bin/ruby

class Phylo

  require 'bio'

  RAXML_BIN = '/home/hoelzer/RubymineProjects/positive_selection/tools/raxml/8.0.25/raxmlHPC-PTHREADS-SSE3'
  NW_DISPLAY_BIN = '/home/hoelzer/RubymineProjects/positive_selection/tools/nw_utilities/nw_display'
  NW_REROOT_BIN = '/home/hoelzer/RubymineProjects/positive_selection/tools/nw_utilities/nw_reroot'
  INKSCAPE_BIN = '/home/hoelzer/RubymineProjects/positive_selection/tools/inkscape'

  attr_accessor :tree_unshod, :tree_corrected, :tree_rooted

  def initialize(alignment, outgroups, out_dir, type, threads, root_species, log, params_string)

    Dir.mkdir(out_dir) unless Dir.exists?(out_dir)

    tree = "#{out_dir}/RAxML_bipartitionsBranchLabels"

    case type
      when :nt
        tree += '.nt'
        if outgroups
          run = "#{RAXML_BIN} -T #{threads} -f a -# 100 -x 1234 -p 1234 -s #{alignment} -n nt -m GTRGAMMA -N 1000 -o #{outgroups} -w #{out_dir}"
          log << run << "\n"
          params_string << run << "\n\n"
          `#{run}`
        else
          run = "#{RAXML_BIN} -T #{threads} -f a -# 100 -x 1234 -p 1234 -s #{alignment} -n nt -m GTRGAMMA -N 1000 -w #{out_dir}"
          log << run << "\n"
          params_string << run << "\n\n"
          `#{run}`
        end
      when :aa
        tree += '.aa'
        if outgroups
          run = "#{RAXML_BIN} -T #{threads} -f a -# 100 -x 1234 -p 1234 -s #{alignment} -n aa -m PROTGAMMAWAG -N 1000 -o #{outgroups} -w #{out_dir}"
          log << run << "\n"
          params_string << run << "\n\n"
          `#{run}`
        else
          run = "#{RAXML_BIN} -T #{threads} -f a -# 100 -x 1234 -p 1234 -s #{alignment} -n aa -m PROTGAMMAWAG -N 1000 -w #{out_dir}"
          log << run << "\n"
          params_string << run << "\n\n"
          `#{run}`
        end
      else
        abort("The chosen type for phylogeny calculation is not valid! #{type}")
    end

    puts "RAXML FINISHED"
    log << "RAXML FINISHED\n\n"

    ##TODO rename the internal species names back to the user defined names to show them in all trees

    # calculate the height of the png tree, for each taxa in the tree reserve 25px
    taxas = 0
    Bio::FastaFormat.open(alignment).each do |e|
      taxas += 1
    end
    export_height = taxas * 25

    @tree_corrected = raxml2drawing(File.open(tree,'r'))

    `#{NW_DISPLAY_BIN} -v 50 -i 'font-size:11' -l 'font-size:16;font-family:helvetica;font-style:italic' -d 'stroke:black;fill:none;stroke-width:2;stroke-linecap:round' -Il -w 650 -b 'opacity:0' -S -s #{tree}.corrected > #{tree}.svg`

    #`#{INKSCAPE_BIN} --verb=FitCanvasToDrawing --verb=FileSave --verb=FileClose #{tree}.svg`
    `#{INKSCAPE_BIN} -f #{tree}.svg -A #{tree}.pdf`

    ## with scale
    cmd = "#{NW_DISPLAY_BIN} -v 50 -i 'font-size:11' -l 'font-size:16;font-family:helvetica;font-style:italic' -d 'stroke:black;fill:none;stroke-width:2;stroke-linecap:round' -Il -w 650 -b 'opacity:0' -s #{tree}.corrected > #{tree}.scale.svg"
    `#{cmd}`
    params_string << cmd << "\n\n"

    #`#{INKSCAPE_BIN} --verb=FitCanvasToDrawing --verb=FileSave --verb=FileClose #{tree}.scale.svg`
    `#{INKSCAPE_BIN} -f #{tree}.scale.svg -A #{tree}.scale.pdf`
    ## png thumbnail
    `#{INKSCAPE_BIN} -f #{tree}.svg --export-width=480 --export-height=#{export_height} --without-gui --export-png=#{tree}.png`

    if root_species
      ## if root species are defined, we can reroot the tree
      cmd = "#{NW_REROOT_BIN} #{tree}.corrected #{root_species.join(' ')} > #{tree}.corrected.rooted"
      `#{cmd}`

      # check if tree was rerooted or LCA error occured
      tree_tmp = File.open("#{tree}.corrected.rooted",'r')
      tree_tmp_string = ''
      tree_tmp.each do |l|
        tree_tmp_string << l
      end
      tree_tmp.close

      if tree_tmp_string.length < 1
        # try reroot with -l option
        cmd = "#{NW_REROOT_BIN} -l #{tree}.corrected #{root_species.join(' ')} > #{tree}.corrected.rooted"
        `#{cmd}`
        params_string << cmd << "\n\n"
        # check again if reroot with -l helped
        tree_tmp = File.open("#{tree}.corrected.rooted",'r')
        tree_tmp_string = ''
        tree_tmp.each do |l|
          tree_tmp_string << l
        end
        tree_tmp.close
      else
        params_string << cmd << "\n\n"
      end

      if tree_tmp_string.length > 1
        `#{NW_DISPLAY_BIN} -v 50 -i 'font-size:11' -l 'font-size:16;font-family:helvetica;font-style:italic' -d 'stroke:black;fill:none;stroke-width:2;stroke-linecap:round' -Il -w 650 -b 'opacity:0' -S -s #{tree}.corrected.rooted > #{tree}.rooted.svg`
        #`#{INKSCAPE_BIN} --verb=FitCanvasToDrawing --verb=FileSave --verb=FileClose #{tree}.rooted.svg`
        `#{INKSCAPE_BIN} -f #{tree}.rooted.svg -A #{tree}.rooted.pdf`
        `#{NW_DISPLAY_BIN} -v 50 -i 'font-size:11' -l 'font-size:16;font-family:helvetica;font-style:italic' -d 'stroke:black;fill:none;stroke-width:2;stroke-linecap:round' -Il -w 650 -b 'opacity:0' -s #{tree}.corrected.rooted > #{tree}.rooted.scale.svg`
        #`#{INKSCAPE_BIN} --verb=FitCanvasToDrawing --verb=FileSave --verb=FileClose #{tree}.rooted.scale.svg`
        `#{INKSCAPE_BIN} -f #{tree}.rooted.scale.svg -A #{tree}.rooted.scale.pdf`
        `#{INKSCAPE_BIN} -f #{tree}.rooted.svg --export-width=480 --export-height=#{export_height} --without-gui --export-png=#{tree}.rooted.png`
      end
      if tree_tmp_string.length > 1
        @tree_rooted = "#{tree}.corrected.rooted"
      else
        @tree_rooted = "#{tree}.corrected"
      end
    end

    # remove bootstrap values from the original tree, they are problematic in CODEML and the file ending has to be .tree!
    # also, use again the internal species names
    tree_file = File.open("#{tree}",'r')
    tree_file_unshod = File.open("#{tree}.unshod.tree",'w')
    tree_file.each do |line|
      tree_file_unshod << line.split(/\[[0-9]+\]/).join('')
    end
    tree_file.close
    tree_file_unshod.close

    @tree_unshod = tree_file_unshod.path

  end

  def raxml2drawing(input)
    output = File.open(input.path+'.corrected','w')

    input.each do |newick|
      # get bootstrap values
      new_newick = newick
      newick.scan(/:[0-9]+\.[0-9]+\[[0-9]+\]/).each do |bootstrap|
        b = bootstrap.split('[')[1].sub(']','')
        n = bootstrap.split('[')[0]
        new_bootstrap = "#{b}#{n}"
        new_newick = new_newick.sub(bootstrap, new_bootstrap)
      end
      # round branch lengths
      newick.scan(/[0-9]+\.[0-9]+/).each do |length|
        f = length.to_f.round(4)
        new_newick = new_newick.sub(length, "#{f}")
      end
      output << new_newick
    end
    output.close
    input.close

    output.path
  end

end