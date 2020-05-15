process tex_built {
    publishDir "${params.output}/${name}/tex/", mode: 'copy', pattern: "*.tex" 
    label 'bioruby'

    input:
        tuple val(name), val(codon_freq), path(mlc), path(aln_aa_nogaps), path(aln_aa_raw), path(internal2input)

    output: 
        tuple val(name), val(codon_freq), path("*.tex"), emit: tex_files
        tuple val(name), path("*_tex", type: 'dir'), emit: tex_dir
        tuple val(name), path("${codon_freq}.lrt"), emit: lrt_params
        tuple val(name), path('gap_start2gap_length.csv'), emit: gap_start2gap_length
        
    script:
    """
    tex.rb ${mlc} ${codon_freq} ${name} ${params.reference} ${aln_aa_nogaps} ${internal2input} ${aln_aa_raw} ${workflow.projectDir}/bin/chi2
    mkdir \$(basename \$PWD)_tex
    cp *.tex \$(basename \$PWD)_tex
    """
}

process frag_tex_built {
    publishDir "${params.output}/${name}/tex/", mode: 'copy', pattern: "*.tex" 
    label 'bioruby'

    input:
        tuple val(name), val(codon_freq), path(mlc), path(aln_aa_nogaps), path(aln_aa_raw), path(internal2input), file(breakpoint_tsv)

    output: 
        tuple val(name), val(codon_freq), path("*.tex"), emit: tex_files
        tuple val(name), path("*_tex", type: 'dir'), emit: tex_dir
        tuple val(name), path("${codon_freq}.lrt"), emit: lrt_params
        
    script:
    """
    tex.rb ${mlc} ${codon_freq} ${name} ${params.reference} ${aln_aa_nogaps} ${internal2input} ${aln_aa_raw} ${workflow.projectDir}/bin/chi2 ${breakpoint_tsv}
    mkdir \$(basename \$PWD)_tex
    cp *.tex \$(basename \$PWD)_tex
    """
}


process pdflatex {
    publishDir "${params.output}/${name}/pdf/", mode: 'copy', pattern: "*.pdf" 
    label 'tex'

    input:
        tuple val(name), val(codon_freq), path(tex_files)

    output: 
        tuple val(name), val(codon_freq), path("*.pdf")
        
    script:
    """
    for TEX in *.tex; do 
        pdflatex \$TEX 
        pdflatex \$TEX # second compilation to adjust the table header width
    done
    """
}

process tex_combine {
    publishDir "${params.output}/${name}/tex/", mode: 'copy', pattern: "*.tex" 
    label 'bioruby'

    input:
        tuple val(name), path(tex_files_1), path(tex_files_2), path(tex_files_3)

    output: 
        tuple val(name), val('full'), path("${name}*.tex")
        
    script:
    """
    tex_combine.rb ${name}_codeml.tex ${name}_gaps_codeml.tex
    """
}
