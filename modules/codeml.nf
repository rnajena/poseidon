process codeml_built {
    label 'bioruby'

    input:
        tuple val(name), path(aln), path(tree)

    output: 
        tuple val(name), path("ctl/*/*/*.ctl")
        
    script:
    """
    codeml_built.rb ctl ${aln} ${tree}
    """
}

process codeml_run {
    label 'codeml'

    input:
        tuple val(name), path(ctl), path(aln), path(tree)

    output: 
        tuple val(name), env(MODEL), path("*.mlc")
        
    script:
    """
    codeml ${ctl} >> ${ctl}.log
    CODON_FREQ=\$(ls ${ctl} | awk 'BEGIN{FS="_"};{print \$2}')
    MODEL=\$(ls ${ctl} | awk 'BEGIN{FS="_"};{print \$3}' | awk 'BEGIN{FS="."};{print \$1}')
    #echo \$CODON_FREQ
    echo \$MODEL
    """
}

// combine all single codeml.mlc files to a final big codeml file
// for each frequency {F61, F1X4, F3X4} build one codeml mlc file

process codeml_combine {
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