process build_fragments {
    label 'bioruby'

    errorStrategy { task.exitStatus = 1 ? 'ignore' :  'terminate' }

    input:
        tuple val(name), path(aln), path(bp), val(recombination)

    output: 
        tuple val(name), path("fragments/fragment_*"), path("fragments/fragment_*/aln/*.nt_ali.checked.nogaps.fasta"), path("fragments/fragment_*/aln/*.aa_ali.checked.nogaps.fasta"), emit: fragments optional true
        tuple val(name), path('bp.tsv'), emit: breakpoints optional true
        
    script:
    """
    if [ ${recombination} == "true" ]; then
        fragments_built.rb ${aln} ${bp}
    fi
    """
}