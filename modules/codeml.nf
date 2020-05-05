process codeml_built {
    publishDir "${params.output}/${name}/codeml/", mode: 'copy', pattern: "ctl/*/*/*.ctl" 
    label 'bioruby'

    input:
        tuple val(name), path(aln), path(tree)

    output: 
        tuple val(name), path("ctl/*/*/*.ctl"), emit: ctl_files
        tuple val(name), path('ctl', type: 'dir'), emit: ctl_dir
        
    script:
    """
    codeml_built.rb ctl ${aln} ${tree}
    """
}

process codeml_run {
    publishDir "${params.output}/${name}/codeml/", mode: 'copy', pattern: "*.mlc" 
    label 'codeml'

    input:
        tuple val(name), path(ctl), path(aln), path(tree)

    output: 
        tuple val(name), env(CODON_FREQ), path("*.mlc"), emit: mlc_files
        
    script:
    """
    codeml ${ctl} >> ${ctl}.log
    CODON_FREQ=\$(ls ${ctl} | awk 'BEGIN{FS="_"};{print \$2}')
    #MODEL=\$(ls ${ctl} | awk 'BEGIN{FS="_"};{print \$3}' | awk 'BEGIN{FS="."};{print \$1}')
    echo \$CODON_FREQ
    """
}

// combine all single codeml.mlc files to a final big codeml file
// for each frequency {F61, F1X4, F3X4} build one codeml mlc file

process codeml_combine {
    publishDir "${params.output}/${name}/codeml/", mode: 'copy', pattern: "*all.mlc" 
    label 'bioruby'

    input:
        tuple val(name), val(codon_freq), path(mlcs)

    output: 
        tuple val(name), val(codon_freq), path("*all.mlc")
        
    script:
    """
    codeml_combine.rb ${codon_freq}
    """
}