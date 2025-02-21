process codeml_built {
    publishDir "${params.output}/${name}/codeml/", mode: 'copy', pattern: "ctl/*/*/*.ctl" 
    label 'bioruby'

    input:
        tuple val(name), path(aln), path(tree)

    output: 
        tuple val(name), path("ctl/*/*/*.ctl"), optional: true, emit: ctl_files
        tuple val(name), path('ctl', type: 'dir'), optional: true, emit: ctl_dir
        
    script:
    """
    codeml_built.rb ctl ${aln} ${tree} ${params.codeml_code}
    #if [ "${aln}" == "x.nt_ali.checked.nogaps.fasta" ]; then rm -r ctl; fi
    """
}

process codeml_run {
    publishDir "${params.output}/${name}/codeml/", mode: 'copy', pattern: "*.mlc" 
    label 'codeml'

    input:
        tuple val(name), path(ctl), path(aln), path(tree)

    output: 
        tuple val(name), env(LABEL), path("*.mlc"), emit: mlc_files
        
    script:
    """
    CODON_FREQ=\$(ls ${ctl} | awk 'BEGIN{FS="_"};{print \$2}')

    #MODEL=\$(ls ${ctl} | awk 'BEGIN{FS="_"};{print \$3}' | awk 'BEGIN{FS="."};{print \$1}')
    #mv *.mlc \$(basename \$PWD).codeml_\${CODON_FREQ}_\${MODEL}.mlc

    LABEL=${name}_\${CODON_FREQ}

    if [ -s ${aln} ]; then
        codeml ${ctl} >> ${ctl}.log
    else
        touch \$(basename \$PWD).dummy.mlc
    fi

    echo \$LABEL
    """
}

// combine all single codeml.mlc files to a final big codeml file
// for each frequency {F61, F1X4, F3X4} build one codeml mlc file

process codeml_combine {
    publishDir "${params.output}/${name}/codeml/", mode: 'copy', pattern: "*all.mlc" 
    label 'bioruby'

    input:
        tuple val(name), val(name_codon_freq), path(mlcs)

    output: 
        tuple val(name), env(CODON_FREQ), path("*all.mlc")
        
    script:
    """
    CODON_FREQ=\$(echo ${name_codon_freq} | sed 's/${name}_//g')
    codeml_combine.rb \${CODON_FREQ}
    """
}

