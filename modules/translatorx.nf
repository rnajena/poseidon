process translatorx {
    label 'translatorx'  
//    publishDir "${params.output}/${name}/", mode: 'copy', pattern: ""

  input:
    tuple val(name), path(fasta)

  output:
    tuple val(name), path("${name}_aln.nt_ali.fasta"), path("${name}_aln.aa_ali.fasta"), path("${name}_aln.aa_based_codon_coloured.html")

  script:
    """
    translatorx -i ${fasta} -p M -o ${name}_aln
    """
  }

process check_aln {
    label 'bioruby'  
//    publishDir "${params.output}/${name}/", mode: 'copy', pattern: ""

  input:
    tuple val(name), path(nt_aln), path(aa_aln), path(html), path(log_file)

  output:
    tuple val(name), path("${name}_aln.nt_ali.checked.fasta"), path("${name}_aln.aa_ali.checked.fasta"), path("${name}_aln.aa_based_codon_coloured.checked.html"), path(log_file)

  script:
    """
    check_aln.rb ${nt_aln} ${aa_aln} ${html} ${name}_aln.nt_ali.checked.fasta ${name}_aln.aa_ali.checked.fasta ${name}_aln.aa_based_codon_coloured.checked.html ${log_file}
    """
  }

process remove_gaps {
    label 'bioruby'  
//    publishDir "${params.output}/${name}/", mode: 'copy', pattern: ""

  input:
    tuple val(name), path(nt_aln), path(aa_aln), path(html), path(log_file)

  output:
    tuple val(name), path("${name}_aln.nt_ali.checked.nogaps.fasta"), emit: fna
    tuple val(name), path("${name}_aln.aa_ali.checked.nogaps.fasta"), emit: faa

  script:
    """
    remove_gaps.rb ${nt_aln} ${aa_aln} ${name}_aln.nt_ali.checked.nogaps.fasta ${name}_aln.aa_ali.checked.nogaps.fasta
    """
  }
