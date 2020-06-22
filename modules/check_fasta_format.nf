process check_fasta_format {
  label 'bioruby'

  input:
    tuple val(name), path(fasta)

  output:
    tuple val(name), path("${name}.checked.fasta"), emit: fasta
    tuple val(name), path("${name}.log"), emit: log_file
    env LOG, emit: log_string 
    tuple val(name), path("${name}_internal2input_species.tsv"), emit: internal2input
    tuple val(name), path("${name}_input2internal_species.tsv"), emit: input2internal
  
  script:
    """
    fasta_format_checker.rb ${fasta} ${name}.checked.fasta ${name}.log ${name}_internal2input_species.tsv ${name}_input2internal_species.tsv ${params.reference} ${params.outgroup}

    # replace newlines with 'foobar' to replace this again in the nextflow for better terminal output
    LOG=\$(sed -z 's/\\n/foobar/g' ${name}.log)
    """
  }
