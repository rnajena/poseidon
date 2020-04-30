process model_selection {
    label 'hyphy'

    input:
        tuple val(name), file(aln), file(trees)

    output: 
        tuple val(name), env(MODEL)
        
    script:
    """
    MODEL='010010' # the default HKY85 model, use this if something happens

    (echo "${aln}"; echo "${name}_nt.raxml.corrected.barefoot.tree"; echo "${params.model_rate_classes}"; echo "${params.model_selection_method}"; echo "${params.model_rejection_level}"; echo "${name}.result") | hyphy mt > model.txt

    HIT=\$(grep 'Model String:' model.txt | awk 'BEGIN{FS=":"};{print \$2}')

    if [ \$HIT ]; then MODEL=\$HIT; fi
    echo \$MODEL
    """
}
