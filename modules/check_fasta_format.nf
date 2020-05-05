process check_fasta_format {
    label 'bioruby'
//    publishDir "${params.output}/${name}/", mode: 'copy', pattern: "dammit"

  input:
    tuple val(name), path(fasta)

  output:
    tuple val(name), path("${name}.checked.fasta"), emit: fasta
    tuple val(name), path("${name}.log"), emit: log
    tuple val(name), path("${name}_internal2input_species.tsv"), emit: internal2input
    tuple val(name), path("${name}_input2internal_species.tsv"), emit: input2internal
  
  script:
    """
    fasta_format_checker.rb ${name} ${fasta} ${name}.checked.fasta ${name}.log ${name}_internal2input_species.tsv ${name}_input2internal_species.tsv ${params.reference} ${params.root}
    """
  }
