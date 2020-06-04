process model_selection {
    label 'hyphy'

    input:
        tuple val(name), file(aln), file(tree)

    output: 
        tuple val(name), env(MODEL), emit: model
        tuple val(name), path('model.log'), emit: log
        
    script:
    """
    MODEL=${params.model} # the default HKY85 model, use this if something happens

    (echo "${aln}"; echo "${tree}"; echo "${params.model_rate_classes}"; echo "${params.model_selection_method}"; echo "${params.model_rejection_level}"; echo "${name}.result") | hyphy mt > model.log

    HIT=\$(grep 'Model String:' model.log | awk 'BEGIN{FS=":"};{print \$2}')

    if [ \$HIT ]; then MODEL=\$HIT; fi
    echo \$MODEL
    """
}

/*
#works
conda install -c bioconda hyphy=2.3.11
cp bats_mx1_* /opt/conda/envs/test/lib/hyphy/TemplateBatchFiles/
(echo "bats_mx1_aln.nt_ali.checked.nogaps.fasta"; echo "bats_mx1_nt.raxml.barefoot.tree"; echo "4"; echo "1"; echo "0.05"; echo "bats_mx1.result") | HYPHYMP /opt/conda/envs/test/lib/hyphy/TemplateBatchFiles/ModelTest.bf

#works
conda install -c bioconda hyphy=2.5.14
(echo "bats_mx1_aln.nt_ali.checked.nogaps.fasta"; echo "bats_mx1_nt.raxml.barefoot.tree"; echo "4"; echo "1"; echo "0.05"; echo "bats_mx1.result") | hyphy mt
*/
