process codeml_built {
    label 'bioruby'

    input:
        tuple val(name), path(aln), path(trees)

    output: 
        tuple val(name), path("ctl/*/*/*.ctl")
        
    script:
    """
    codeml_built.rb ctl ${aln} ${name}_nt.raxml.corrected.barefoot.tree
    """
}

process codeml_run {
    label 'codeml'

    input:
        tuple val(name), path(ctl), path(aln), path(trees)

//    output: 
//        tuple val(name), env(MODEL)
        
    script:
    """
    codeml ${ctl} >> ${ctl}.log
    """
}
