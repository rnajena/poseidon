process gard_detect {
    label 'gard'

    input:
        tuple val(name), file(aln), val(model)

    output: 
        tuple val(name), path(aln), path("gard.html"), path("gard.log"), path("gard_processor.log")
        
    script:
    """
    GARD_TEMPLATE_BATCH='/usr/lib/hyphy/TemplateBatchFiles/GARD.bf'
    GARD_PROCESSOR_TEMPLATE_BATCH='/usr/lib/hyphy/TemplateBatchFiles/GARDProcessor.bf'
    OPENMPI='mpirun'
    OPENMPI_RUN='--allow-run-as-root'
    HYPHYMPI='hyphympi'

    TMPDIR=\${PWD}

    OUTPUT=\${TMPDIR}/gard
    ALN=\${TMPDIR}/${aln}
    MODEL=${model}

    (echo "\${ALN}"; echo "\${ALN}"; echo "\${MODEL}"; echo "${params.gard_rate_variation}"; echo "${params.gard_rate_classes}"; echo "\${OUTPUT}") | \${OPENMPI} \${OPENMPI_RUN} -np ${task.cpus} \${HYPHYMPI} \${GARD_TEMPLATE_BATCH} &> \${TMPDIR}/gard.log

    mv gard gard.html

    # check if GARD was able to detect breakpoints
    if grep -q "ERROR: Too few sites for c-AIC inference." gard.log; then 
        touch gard_processor.log
        echo 'Warning: PoSeiDon was not able to perform a recombination analysis with GARD because of too few sites for c-AIC inference. We do not recommend to use this alignment for positive selection detection. Please be careful with any sites detected as positively selected! If possible, extend your sequences or reduce the number of entries in your FASTA file.' > gard_processor.log
    else
        if grep -q "ERROR: There are too few potential break points to support" gard.html; then
            touch gard_processor.log
            echo 'Warning: PoSeiDon was not able to perform a recombination analysis with GARD because there are too few potential break points to support recombination events. We do not recommend to use this alignment for positive selection detection. Please be careful with any sites detected as positively selected! If possible, extend your sequences or the number of entries in your FASTA file.' > gard_processor.log
        else 
            echo 'Found breakpoints! Run processing and KH test!'

            gard_result_file=\${TMPDIR}/gard_finalout
            gard_splits_file=\${TMPDIR}/gard_splits
            gard_processor_log=\${TMPDIR}/gard_processor.log
            (echo "\${gard_result_file}"; echo "\${gard_splits_file}") | \${OPENMPI} \${OPENMPI_RUN} -np ${task.cpus} \${HYPHYMPI} \${GARD_PROCESSOR_TEMPLATE_BATCH} &> \$gard_processor_log
        fi
    fi

    """
}

process gard_process {
    label 'bioruby'
    publishDir "${params.output}/${name}/recombination/", mode: 'copy', pattern: "gard.adjusted.html"
    
    input:
        tuple val(name), path(aln), path(gard_html), path(gard_log), path(gard_processor_log)

    output: 
        tuple val(name), path("gard.adjusted.html"), emit: html
        tuple val(name), path("bp.tsv"), emit: bp
        tuple val(name), env(RECOMBINATION), emit: recombination
        
    script:
    """
    # check if GARD was able to detect breakpoints
    if grep -q "Warning: PoSeiDon was not able to perform a recombination analysis with GARD" gard_processor.log; then 
        # nothing, maybe write dummy files for output channel
        touch gard.adjusted.html
        touch bp.tsv
        cat gard_processor.log > gard.adjusted.html 
        RECOMBINATION=false
    else     
        RECOMBINATION=true
        gard.rb ${params.kh} ${gard_html} ${aln}
        if [[ \$(grep -v '#' bp.tsv | wc -l) == 0 ]]; then
            RECOMBINATION=false
        fi
    fi
    echo \$RECOMBINATION
    """
}
