process raxml_nt {
    label 'raxml'  
//    publishDir "${params.output}/${name}/", mode: 'copy', pattern: ""

  input: 
    tuple val(name), file(aln)

  output:
    tuple val(name), file("${name}_nt.raxml")

  script:
    if (params.root == 'NA')
    """
      raxmlHPC-PTHREADS-SSE3 -T ${task.cpus} -f a -x 1234 -p 1234 -s ${aln} -n nt -m GTRGAMMA -N ${params.bootstrap} 
      mv RAxML_bipartitionsBranchLabels.nt ${name}_nt.raxml
    """
    else
    """
      raxmlHPC-PTHREADS-SSE3 -T ${task.cpus} -f a -x 1234 -p 1234 -s ${aln} -n nt -m GTRGAMMA -N ${params.bootstrap} -o ${params.root}
      mv RAxML_bipartitionsBranchLabels.nt ${name}_nt.raxml
    """    
}

process raxml_aa {
    label 'raxml'  
//    publishDir "${params.output}/${name}/", mode: 'copy', pattern: ""

  input: 
    tuple val(name), file(aln)

  output:
    tuple val(name), file("${name}_aa.raxml")

  script:
    if (params.root == 'NA')
    """
      raxmlHPC-PTHREADS-SSE3 -T ${task.cpus} -f a -x 1234 -p 1234 -s ${aln} -n aa -m PROTGAMMAWAG -N ${params.bootstrap} 
      mv RAxML_bipartitionsBranchLabels.aa ${name}_aa.raxml
    """
    else
    """
      raxmlHPC-PTHREADS-SSE3 -T ${task.cpus} -f a -x 1234 -p 1234 -s ${aln} -n aa -m PROTGAMMAWAG -N ${params.bootstrap} -o ${params.root}
      mv RAxML_bipartitionsBranchLabels.aa ${name}_aa.raxml
    """    
}

process raxml2drawing {
  label 'bioruby'

  input:
    tuple val(name), file(raxml_nt), file(raxml_aa)

  output:
    tuple val(name), file("${raxml_nt.baseName}.raxml.corrected"), file("${raxml_aa.baseName}.raxml.corrected")

  script:
  """
    raxml2drawing.rb ${raxml_nt} ${raxml_nt.baseName}.raxml.corrected
    raxml2drawing.rb ${raxml_aa} ${raxml_aa.baseName}.raxml.corrected
  """
}

/*
TODO: it seems that only one inkscape command can be started?!?! test maxForks 1
*/
//bats_mx1_nt.raxml.corrected.pdf  bats_mx1_nt.raxml.corrected.scale.pdf  bats_mx1_nt.raxml.corrected.svg
//bats_mx1_nt.raxml.corrected.png  bats_mx1_nt.raxml.corrected.scale.svg
process nw_display {
  label 'newick_utils'

  maxForks 1

  input:
    tuple val(name), file(newick_nt), file(newick_aa), file(aln)

  output:
    tuple val(name), file("*.svg"), file("*.pdf"), file("*.png")

  script:
  """
    for newick in *.corrected; do 

    nw_display -v 50 -i 'font-size:11' -l 'font-size:16;font-family:helvetica;font-style:italic' -d 'stroke:black;fill:none;stroke-width:2;stroke-linecap:round' -Il -w 650 -b 'opacity:0' -S -s \${newick} > \${newick}.svg
    inkscape -f \${newick}.svg -A \${newick}.pdf 
    #convert \${newick}.svg \${newick}.pdf

    nw_display -v 50 -i 'font-size:11' -l 'font-size:16;font-family:helvetica;font-style:italic' -d 'stroke:black;fill:none;stroke-width:2;stroke-linecap:round' -Il -w 650 -b 'opacity:0' -s \${newick} > \${newick}.scale.svg
    inkscape -f \${newick}.scale.svg -A \${newick}.scale.pdf

    # png thumbnail
    ## calculate the export height based on the number of taxa
    TAXA=\$(grep ">" ${aln} | wc -l | awk '{print \$1}')
    HEIGHT=\$((25*TAXA))
    inkscape -f \${newick}.svg --export-width=480 --export-height=\$HEIGHT --without-gui --export-png=\${newick}.png

    done
  """
}

/*
remove bootstrap values from the tree, they are problematic in CODEML and the file ending has to be .tree!
*/
process barefoot {
  label 'bioruby'

  input:
    tuple val(name), file(newick)

  output:
    tuple val(name), file("*.barefoot.tree")

  script:
    """
    #!/usr/bin/env ruby

    tree_file = File.open("${newick}",'r')
    tree_file_barefoot = File.open("${newick}.barefoot.tree",'w')
    tree_file.each do |line|
      tree_file_barefoot << line.split(/\\[[0-9]+\\]/).join('')
    end
    tree_file.close
    tree_file_barefoot.close
    """
}

/*
TODO: the second part (the drawing) can be done with the already available function
from phylo.rb
// IMPORTANT: file should be named with .rooted then: bats_mx1_nt.raxml.corrected.rooted for further processing
*/
/*process reroot {
  
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

 } */


