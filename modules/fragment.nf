process build_fragments {
    label 'bioruby'

    errorStrategy { task.exitStatus = 1 ? 'ignore' :  'terminate' }

    input:
        tuple val(name), path(aln), path(bp), val(recombination)

    output: 
        tuple val(name), path("fragments/fragment_*"), emit: fragments
        
    script:
    """
    if [ ${recombination} == "true" ]; then
        fragments_built.rb ${aln} ${bp}
    else
        echo false
    fi
    """
}
