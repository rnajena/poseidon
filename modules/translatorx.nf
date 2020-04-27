process dammit {
    label 'dammit'  
    publishDir "${params.output}/${name}/", mode: 'copy', pattern: "dammit"

  input:
    tuple val(name), path(transcriptome)
    path(dbs)

  output:
    tuple val(name), path("dammit")

  script:
    if (params.full)
    """
    BUSCO=\$(echo ${params.busco} | awk 'BEGIN{FS="_"};{print \$1}')
    dammit annotate ${transcriptome} --database-dir \${PWD}/dbs --busco-group \${BUSCO} -n ${name} -o dammit --n_threads ${task.cpus} --full 
    """
    else
    """
    BUSCO=\$(echo ${params.busco} | awk 'BEGIN{FS="_"};{print \$1}')
    dammit annotate ${transcriptome} --database-dir \${PWD}/dbs --busco-group \${BUSCO} -n ${name} -o dammit --n_threads ${task.cpus} #--full 
    """
  }
