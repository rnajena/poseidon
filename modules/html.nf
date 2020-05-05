process html_main {
    //publishDir "${params.output}/${name}/tex/", mode: 'copy', pattern: "*.tex" 
    label 'bioruby'

    input:
        val type 
        tuple val(name), path(aa_aln_html), path(aa_aln), path(ctl_dir), path(tree_nt), path(tree_aa), path(input_fasta), path(internal2input_species), path(input2internal_species), path(gard_html_file), val(nucleotide_bias_model), path(tex_files), path(tex_dir), path(mlc_files_1), path(mlc_files_2), path(mlc_files_3), path(lrt_files), path(pdfs), path(nt_aln_checked)

//    output: 
//        tuple val(name), val(codon_freq), path("*.tex")
        
    script:
    """
    mkdir -p html/full_aln
    touch html/full_aln/index.html
    ADJUSTED_DOMAIN_POS='NA'
    TITLE=${name}
    ALN_LENGTH_WITH_GAPS='0'

    # copy mlc files in their correct ctl folder
    # codeml_F61_M8a.mlc
    # No such file or directory @ rb_sysopen - ctl/F3X4/M0/codeml.mlc
    mkdir ctl_mlc
    for MLC in *.mlc; do
        FREQ=\$(basename \$MLC .mlc | awk 'BEGIN{FS="_"};{print \$2}')
        MODEL=\$(basename \$MLC .mlc | awk 'BEGIN{FS="_"};{print \$3}')
        mkdir -p ctl_mlc/\$FREQ/\$MODEL
        cp \$MLC ctl_mlc/\$FREQ/\$MODEL/codeml.mlc
        cp ctl/\$FREQ/\$MODEL/*.ctl ctl_mlc/\$FREQ/\$MODEL/codeml.clt
    done

    mkdir tex
    cp *_tex/*.tex tex/

    html.rb ${type} html html/full_aln/index.html ${aa_aln_html} ${aa_aln} ctl_mlc ${tree_nt} ${tree_aa} \$ADJUSTED_DOMAIN_POS \$TITLE ${input_fasta} ${internal2input_species} ${input2internal_species} \$ALN_LENGTH_WITH_GAPS ${gard_html_file} ${nucleotide_bias_model} ${name}_codeml.tex ${name}_gaps_codeml.tex tex ${params.refactor} ${params.poseidon_version} ${params.recombination} ${workflow.projectDir}
    """
}

// TODO: add the mlc files to the CTL dir!
// TODO: check params.refactor
// TODO: check params.recombination
// TODO: FRAGMENTS!!

/*
type, html_dir, out, translatorx_html, aa_aln, codeml_results, nt_tree, aa_tree, domain_pos, title, input_fasta, internal2input_species, input2internal_species, aln_length_with_gaps_adjustor, gard_html_file, nucleotide_bias_model, index_html_paths, tex_summary_file_path, tex_summary_file_path_gapped, tex_objects, refactored_aln, version, is_recomb
*/
