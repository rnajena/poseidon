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
        tuple val(name), path('gap_start2gap_length.csv'), emit: gap_start2gap_length
        
    script:
    """
    if [ -s ${mlc} ]; then
        tex.rb ${mlc} ${codon_freq} ${name} ${params.reference} ${aln_aa_nogaps} ${internal2input} ${aln_aa_raw} ${workflow.projectDir}/bin/chi2 ${breakpoint_tsv}
        mkdir \$(basename \$PWD)_tex
        cp *.tex \$(basename \$PWD)_tex
    else
        touch dummy.tex
        mkdir \$(basename \$PWD)_dummy_tex
        touch ${codon_freq}.lrt
        touch gap_start2gap_length.csv
    fi
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
    if [ ! -f dummy.tex ]; then
        for TEX in *.tex; do 
            pdflatex \$TEX 
            pdflatex \$TEX # second compilation to adjust the table header width
        done
    else
        touch dummy.pdf
    fi
    """
}

process tex_combine {
    publishDir "${params.output}/${name}/tex/", mode: 'copy', pattern: "*.tex" 
    label 'bioruby'

    input:
//        tuple val(name), path(tex_dir_1), path(tex_dir_2), path(tex_dir_3)
        tuple val(name), path(full_aln_tex_dirs), path(frag_tex_dirs), file(bp_tsv)

    output: 
        tuple val(name), val('full'), path("${name}*.tex")
        
    script:
    """
    cp *_tex/* .
    tex_combine.rb ${name}_codeml.tex ${name}_gaps_codeml.tex ${bp_tsv}
    """
}
