process model_selection {
    label 'hyphy'

    input:
        tuple val(name), file(aln), file(tree)

    output: 
        tuple val(name), env(MODEL), emit: model
        tuple val(name), path('model.log'), emit: log
        
    script:
    """
    MODEL='010010' # the default HKY85 model, use this if something happens

    (echo "${aln}"; echo "${tree}"; echo "${params.model_rate_classes}"; echo "${params.model_selection_method}"; echo "${params.model_rejection_level}"; echo "${name}.result") | hyphy mt > model.log

    HIT=\$(grep 'Model String:' model.log | awk 'BEGIN{FS=":"};{print \$2}')

    if [ \$HIT ]; then MODEL=\$HIT; fi
    echo \$MODEL
    """
}
