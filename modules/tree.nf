process raxml_nt {
    label 'raxml'  
//    publishDir "${params.output}/${name}/", mode: 'copy', pattern: ""

  errorStrategy { task.exitStatus = 255 ? 'ignore' : 'terminate' }

  input: 
    tuple val(name), file(aln)

  output:
    tuple val(name), file("${name}_nt.raxml")

  script:
    if (params.outgroup == 'NA')
    """
      raxmlHPC-PTHREADS-SSE3 -T ${task.cpus} -f a -x 1234 -p 1234 -s ${aln} -n nt -m GTRGAMMA -N ${params.bootstrap} 
      mv RAxML_bipartitionsBranchLabels.nt ${name}_nt.raxml
    """
    else
    """
      raxmlHPC-PTHREADS-SSE3 -T ${task.cpus} -f a -x 1234 -p 1234 -s ${aln} -n nt -m GTRGAMMA -N ${params.bootstrap} -o ${params.outgroup.toUpperCase()}
      mv RAxML_bipartitionsBranchLabels.nt ${name}_nt.raxml
    """    
}

process raxml_aa {
    label 'raxml'  
//    publishDir "${params.output}/${name}/", mode: 'copy', pattern: ""

  errorStrategy { task.exitStatus = 255 ? 'ignore' : 'terminate' }

  input: 
    tuple val(name), file(aln)

  output:
    tuple val(name), file("${name}_aa.raxml")

  script:
    if (params.outgroup == 'NA')
    """
      raxmlHPC-PTHREADS-SSE3 -T ${task.cpus} -f a -x 1234 -p 1234 -s ${aln} -n aa -m PROTGAMMAWAG -N ${params.bootstrap} 
      mv RAxML_bipartitionsBranchLabels.aa ${name}_aa.raxml
    """
    else
    """
      raxmlHPC-PTHREADS-SSE3 -T ${task.cpus} -f a -x 1234 -p 1234 -s ${aln} -n aa -m PROTGAMMAWAG -N ${params.bootstrap} -o ${params.outgroup.toUpperCase()}
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
    for newick in *.corrected*; do 

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
process reroot {
  label 'newick_utils'

  input:
    tuple val(name), path(nt_newick), path(aa_newick)

  output:
    tuple val(name), file("${nt_newick}.rooted"), emit: nt optional true
    tuple val(name), file("${aa_newick}.rooted"), emit: aa optional true
    tuple val(name), env(ROOTED), emit: worked

  // TODO: a comma separated list can be in params.outgroup! split and join(' ')!
  script:  
    """
    ROOTED=false

    # if root species are defined, we can reroot the tree
    nw_reroot ${nt_newick} ${params.outgroup.toUpperCase()} > ${nt_newick}.rooted
    nw_reroot ${aa_newick} ${params.outgroup.toUpperCase()} > ${aa_newick}.rooted

    # check if tree was rerooted or LCA error occured
    LINES=\$(wc -l ${nt_newick}.rooted | awk '{print \$1}')

    if [[ \$LINES < 1 ]]; then 
        # try reroot with -l option
        nw_reroot -l ${nt_newick} ${params.outgroup.toUpperCase()} > ${nt_newick}.rooted
        nw_reroot -l ${aa_newick} ${params.outgroup.toUpperCase()} > ${aa_newick}.rooted
        LINES=\$(wc -l ${nt_newick}.rooted | awk '{print \$1}')
    fi

    if [[ \$LINES > 1 ]]; then 
      ROOTED=true
    fi
    """
 }


