process tex_built {
    label 'bioruby'

    input:
        tuple val(name), val(codon_freq), path(mlc), path(aln_aa_nogaps), path(aln_aa_raw), path(internal2input)

    output: 
        tuple val(name), val(codon_freq), path("*.tex")
        
    script:
    """
    tex.rb ${mlc} ${codon_freq} ${name} ${params.reference} ${aln_aa_nogaps} ${internal2input} ${aln_aa_raw} ${workflow.projectDir}/bin/chi2
    """
}

process pdflatex {
    label 'tex'

    input:
        tuple val(name), val(codon_freq), path(tex_files)

    output: 
        tuple val(name), val(codon_freq), path("*.pdf")
        
    script:
    """
    for TEX in *.tex; do pdflatex \$TEX; done
    """
}
