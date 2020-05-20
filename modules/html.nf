process html {
    //if (type == 'full_aln')
        publishDir "${params.output}/${name}/", mode: 'copy', pattern: "html" 
    //else
        publishDir "${params.output}/${name.split('_fragment_')[0]}/html/", mode: 'copy', pattern: "fragment_*" 

    label 'bioruby'

    maxForks 1

    input:
        val type 
        tuple val(name), path(aa_aln_html), path(aa_aln), path(ctl_dir), path(tree_nt), path(tree_aa), path(input_fasta), path(internal2input_species), path(input2internal_species), path(gard_html_file), val(nucleotide_bias_model), path(tex_files), path(tex_dir), path(mlc_files_1), path(mlc_files_2), path(mlc_files_3), path(lrt_files), path(pdfs), path(nt_aln_checked), path(tree_svg), path(tree_pdf), path(tree_png), path(model_log), path(translated_fasta), path(aln_nt_nogaps), path(aln_aa_nogaps), file(aln_length_with_gaps), file(gap_adjusted_start_end), val(recombination)

    output: 
        tuple val(name), file("html/*/index.html"), emit: index
        tuple val(name), file("tex"), emit: tex_dir
        path("html", type: 'dir')
        path("fragment_*", type: 'dir') optional true
        
    script:
    """
    if [[ ${type} == 'fragment' ]]; then
        # extract the fragment ID from the name_frag combination label
        NAME=\$(echo ${name} | awk 'BEGIN{FS="_fragment_"};{print "fragment_"\$2}')
        FASTA_NAME=\$(echo ${name} | awk 'BEGIN{FS="_fragment_"};{print \$1}')
        TYPE=\$NAME
    else
        NAME=${name}
        FASTA_NAME=${name}
        TYPE='full_aln'
    fi

    mkdir -p html/\${TYPE}
    touch html/\${TYPE}/index.html
    TITLE=\${NAME}

    # copy mlc files in their correct ctl folder
    # codeml_F61_M8a.mlc
    # No such file or directory @ rb_sysopen - ctl/F3X4/M0/codeml.mlc
    mkdir ctl_mlc
    for MLC in *.mlc; do
        FREQ=\$(basename \$MLC .mlc | awk 'BEGIN{FS="_"};{print \$2}')
        MODEL=\$(basename \$MLC .mlc | awk 'BEGIN{FS="_"};{print \$3}')
        mkdir -p ctl_mlc/\$FREQ/\$MODEL
        cp \$MLC ctl_mlc/\$FREQ/\$MODEL/codeml.mlc
        cp ctl/\$FREQ/\$MODEL/*.ctl ctl_mlc/\$FREQ/\$MODEL/codeml.ctl
    done

    mkdir tex
    cp *_tex/*.tex tex/

    html.rb \${TYPE} html html/\${TYPE}/index.html ${aa_aln_html} ${aa_aln} ctl_mlc ${tree_nt} ${tree_aa} ${gap_adjusted_start_end} \$TITLE ${input_fasta} ${internal2input_species} ${input2internal_species} ${aln_length_with_gaps} ${gard_html_file} ${nucleotide_bias_model} \${FASTA_NAME}_codeml.tex \${FASTA_NAME}_gaps_codeml.tex tex ${params.refactor} ${params.poseidon_version} ${recombination} ${workflow.projectDir}

    if [[ ${type} == 'fragment' ]]; then
        cp -r html/\$NAME .
    fi

    sleep 5s
    """
}


// TODO: add the mlc files to the CTL dir!
// TODO: check params.refactor

/*
type, html_dir, out, translatorx_html, aa_aln, codeml_results, nt_tree, aa_tree, domain_pos, title, input_fasta, internal2input_species, input2internal_species, aln_length_with_gaps_adjustor, gard_html_file, nucleotide_bias_model, index_html_paths, tex_summary_file_path, tex_summary_file_path_gapped, tex_objects, refactored_aln, version, is_recomb
*/



process html_codeml {
    publishDir "${params.output}/${name}/html/${type}/", mode: 'copy', pattern: "codeml.html" 
    label 'bioruby'

    input:
        val type 
        tuple val(name), file(html_main_index), path(tex_dir), path(lrt_files), val(fragment_names)

    output: 
        tuple val(name), path("codeml.html")
        
    script:
    """
    codeml_html.rb ${type} codeml.html ${html_main_index} ${tex_dir} '${fragment_names}'
    """
}



/*parameter_html_out, html_index_file, frag_names, timestamp, version, parameter_strings*/
process html_params {
    publishDir "${params.output}/${name}/html/${type}/", mode: 'copy', pattern: "params.html" 
    label 'bioruby'

    input:
        val type 
        tuple val(name), path(html_main_index)

    output: 
        tuple val(name), path("params.html")
        
    script:
    """
    parameter_html.rb params.html ${html_main_index} ${params.poseidon_version} ${workflow.projectDir}
    """
}


/**/
process frag_aln_html {
//    publishDir "${params.output}/${name}/html/${type}/", mode: 'copy', pattern: "params.html" 
    label 'bioruby'

    input:
        tuple val(name_frag), val(name), val(frag), path(gap_start2gap_length_csv), path(aa_aln_html_raw), path(bp_tsv), path(nt_aln_raw), path(aa_aln_raw)

    output: 
        tuple val(name_frag), val(name), val(frag), path("fragments/${frag}/aln/${name}_aln.aa_based_codon_coloured.html"), path("fragments/${frag}/aln/${name}_aln.aa_ali.fasta"), path("fragments/${frag}/aln/${name}_aln.nt_ali.fasta"), emit: all 
        tuple val(name_frag), path("${frag}_gap-adjusted_start_end.csv"), emit: gap_adjusted_start_end
        tuple val(name_frag), path('aa_bp_with_gaps.csv'), emit: aa_bp_with_gaps
        //tuple env(NEXT_NAME_FRAG), env(ALN_LENGTH_WITH_GAPS), emit: aln_length_with_gaps
        //tuple env(DUMMY_NAME_FRAG), val('1'), emit: dummy // for the first fragment
    script:
    """
    mkdir -p fragments/${frag}/aln
    ALN_LENGTH_WITH_GAPS=\$(frag_aln_html.rb ${aa_aln_html_raw} ${frag} ${bp_tsv} ${aa_aln_raw} ${gap_start2gap_length_csv} | awk 'BEGIN{FS="aln_length_with_gaps:"}{print \$2}')

    # the calculated aln_length_with_gaps are needed for the NEXT fragment
    FRAG_COUNT=\$(echo ${frag} | sed 's/fragment_//g')
    FRAG_COUNT=\$((\$FRAG_COUNT+1))
    NEXT_FRAG="fragment_\${FRAG_COUNT}"
    NEXT_NAME_FRAG="${name}_\${NEXT_FRAG}" # bats_mx1_fragment_1 -> bats_mx1_fragment_2
    DUMMY_NAME_FRAG="${name}_fragment_1"

    echo \$ALN_LENGTH_WITH_GAPS
    """
}


process html_recomb {
    publishDir "${params.output}/${name}/html/full_aln/", mode: 'copy', pattern: "recomb.html" 
    label 'bioruby'

    input:
        tuple val(name), path(html_index), path(gard_html), path(full_nt_tree), path(frag_nt_trees)

    output: 
        path "recomb.html"
        
    script:
    """
    mkdir -p html/full_aln/
    touch html/full_aln/recomb.html
    recombination_html.rb ${html_index} ${gard_html} ${full_nt_tree} '${frag_nt_trees}'
    cp html/full_aln/recomb.html recomb.html
    """
}