process gard {
    label 'gard'

    input:
        tuple val(name), file(aln), val(model)

//    output: 
//        tuple val(name), env(MODEL)
        
    script:
    """
    # it seems that a global varibale is simply used by HYPHY/GARD
    CPU=${task.cpus}

    (echo "${aln}"; echo "${aln}"; echo "${model}"; echo "${params.gard_rate_variation}"; echo "${params.gard_rate_classes}"; echo "gard.out") | hyphy gard > gard.txt

    """
}
