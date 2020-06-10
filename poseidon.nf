#!/usr/bin/env nextflow
nextflow.preview.dsl=2

/*
* PoSeiDon -- Positive Selection Detection and Recombination Analysis
* Author: hoelzer.martin@gmail.com
*/

/************************** 
* Help messages, user inputs & checks
**************************/

if( !nextflow.version.matches('20.+') ) {
    println "This workflow requires Nextflow version 20.X or greater -- You are running version $nextflow.version"
    exit 1
}

if (params.help) { exit 0, helpMSG() }

println " "
println "\u001B[32mProfile: $workflow.profile\033[0m"
println " "
println "\033[2mCurrent User: $workflow.userName"
println "Nextflow-version: $nextflow.version"
println "Starting time: $nextflow.timestamp"
println "Workdir location:"
println "  $workflow.workDir\u001B[0m"
println " "
if (workflow.profile == 'standard') {
println "\033[2mCPUs to use: $params.cores"
println "Output dir name: $params.output\u001B[0m"
println " "}

println "\033[2mTree root species: $params.outgroup"
println "Reference species: $params.reference\u001B[0m"
println " "
if (params.kh) {
    println "\033[2mUse KH-insignificant breakpoints: yes\u001B[0m"
    println " "
} else {
    println "\033[2mUse KH-insignificant breakpoints: no\u001B[0m"
    println " "
}

if (params.profile) {
    exit 1, "--profile is WRONG use -profile" }
if (!params.fasta) {
    exit 1, "input missing, use [--fasta]"}
    
// fasta input
    if (params.fasta && params.list) { fasta_input_ch = Channel
            .fromPath( params.fasta, checkIfExists: true )
            .splitCsv()
            .map { row -> [row[0], file("${row[1]}", checkIfExists: true)] }
            .view() }
    else if (params.fasta) { fasta_input_ch = Channel
            .fromPath( params.fasta, checkIfExists: true, type: 'file')
            .map { file -> tuple(file.simpleName, file) }
            .view() }

/************************** 
* MODULES
**************************/

include check_fasta_format from './modules/check_fasta_format'
include {translatorx; check_aln; remove_gaps} from './modules/translatorx'

include {raxml_nt; raxml_aa; raxml2drawing; nw_display; barefoot; reroot; reroot as frag_reroot} from './modules/tree'
include raxml_nt as frag_raxml_nt from './modules/tree'
include raxml_aa as frag_raxml_aa from './modules/tree'
include raxml2drawing as frag_raxml2drawing from './modules/tree'
include nw_display as frag_nw_display from './modules/tree'
include barefoot as frag_barefoot from './modules/tree'

include model_selection from './modules/model_selection'

include {gard_detect; gard_process} from './modules/gard'

include {codeml_run; codeml_built; codeml_combine} from './modules/codeml'
include codeml_run as frag_codeml_run from './modules/codeml'
include codeml_built as frag_codeml_built from './modules/codeml'
include codeml_combine as frag_codeml_combine from './modules/codeml'

include {tex_built; tex_combine} from './modules/tex'
include frag_tex_built from './modules/tex'
include pdflatex as pdflatex_single from './modules/tex'
include pdflatex as frag_pdflatex_single from './modules/tex'
include pdflatex as pdflatex_full from './modules/tex'

include {html; html as html_frag; html_codeml; html_codeml as frag_html_codeml; html_params; frag_aln_html; html_recomb} from './modules/html'

include build_fragments from './modules/fragment'


/************************** 
* MAIN WORKFLOW 
**************************/

workflow {

    /*******************************
    check input FASTA */
    check_fasta_format(fasta_input_ch)
    
    /*******************************
    TODO rndm subsampling if too many sequences, PAML is not inteded for >100 sequences */

    /*******************************
    align, check alignment, and remove gaps */
    remove_gaps(
        check_aln(
            translatorx(
                check_fasta_format.out.fasta).all
                .join(check_fasta_format.out.log)
        )
    )

    /*******************************
    raw alignments and html from translatorx */
    raw_aln_aa = translatorx.out.all.map {name, nt_aln, aa_aln, html -> [name, aa_aln] }
    raw_aln_nt = translatorx.out.all.map {name, nt_aln, aa_aln, html -> [name, nt_aln] }
    raw_aln_html = translatorx.out.all.map {name, nt_aln, aa_aln, html -> [name, html] }

    /*******************************
    checked alignments and html from translatorx */
    checked_aln_aa = check_aln.out.map {name, nt_aln, aa_aln, html, log -> [name, aa_aln] }
    checked_aln_nt = check_aln.out.map {name, nt_aln, aa_aln, html, log -> [name, nt_aln] }
    checked_aln_html = check_aln.out.map {name, nt_aln, aa_aln, html, log -> [name, html] }

    /*******************************
    ID mappings */ 
    internal2input_c = check_fasta_format.out.internal2input
    input2internal_c = check_fasta_format.out.input2internal

    /*******************************
    build nt and aa phylogeny */
    tree_nt = raxml_nt(remove_gaps.out.fna)
    tree_aa = raxml_aa(remove_gaps.out.faa)

    /*******************************
    convert RAxML output newicks for drawing */
    newicks_ch = 
        raxml2drawing(
            raxml_nt.out.combine(raxml_aa.out, by: 0)
        ).combine(
            remove_gaps.out.fna, by: 0 // simply needed to count the number of taxa for plotting height estimation
        )

    /*******************************
    if outgroups are provided, try reroot the trees */
    if (params.outgroup != 'NA' ) {
        reroot(newicks_ch.map{name, nt_tree, aa_tree, nt_nogaps_fasta -> tuple (name, nt_tree, aa_tree)})
        tree_nt = reroot.out.nt
        tree_aa = reroot.out.aa
        // update the newicks_ch with the re-rooted trees
        newicks_ch = tree_nt.join(tree_aa, by: 0).combine(remove_gaps.out.fna, by: 0)
    }

    /*******************************
    draw trees */
    nw_display(newicks_ch)

    /*******************************
    full nt aln w/o gaps and the corresponding nt tree w/o bootstraps */
    aln_tree_ch = remove_gaps.out.fna.join(barefoot(raxml_nt.out))

    /*******************************
    model selection */
    model_selection(aln_tree_ch)

    /*******************************
    GARD breakpoint analysis */
    gard_process(
        gard_detect(remove_gaps.out.fna.join(model_selection.out.model))
    )

    /*******************************
    CODEML */
    codeml_combine(
            codeml_run(
                codeml_built(aln_tree_ch).ctl_files.transpose().combine(aln_tree_ch, by: 0)
            ).mlc_files.groupTuple(by: 1, size: 6).map { it -> tuple ( it[0][1], it[1], it[2])} //[bats_mx1, bats_mx1_F1X4, [codeml_F1X4_M1a.mlc, codeml_F1X4_M0.mlc, codeml_F1X4_M8a.mlc, codeml_F1X4_M7.mlc, codeml_F1X4_M8.mlc, codeml_F1X4_M2a.mlc]]
    ) //[bats_mx1, F61, /home/martin/git/poseidon/work/7b/43745a2f2434cd51c0ea91945261e6/codeml_F61.all.mlc]

    /*******************************
    LaTeX summary tables for single comparisons */
    pdflatex_single(
            tex_built(
                codeml_combine.out
                    .combine(remove_gaps.out.faa, by: 0)
                    .combine(raw_aln_aa, by: 0)
                    .combine(internal2input_c, by: 0)
            ).tex_files
    )

    // this is needed to plot the color bars if there are fragments
    // channel will be filled in fragment part if there are fragments and then used in full html module
    // otherwise, this dummy channel will tell the full html module to not plot any color bars, 'NA'
    // if there is a fragment file providing information about fragment positions, this will be used for the full_aln
    // in the html.rb script
    gap_adjusted_start_end_main = checked_aln_html.map{name, html -> tuple (name, 'NA')}
    frag_aa_bp_with_gaps_for_full_aln = checked_aln_html.map{name, html -> tuple (name, '1')}
    fragment_names_c = checked_aln_html.map{name, html -> tuple (name, 'NA')}

    /*******************************
    -- FRAGMENTS --
    If breakpoints ... */ 
    fragments_ch = build_fragments(remove_gaps.out.fna.join(gard_process.out.bp).join(gard_process.out.recombination)).fragments.transpose()

    // for each fragment, build a new unrooted tree for CODEML
    frag_aln_nogaps_nt = fragments_ch.map{ name, fragment_dir, fragment_nt, fragment_aa -> tuple ("${name}_${fragment_dir.getFileName()}", fragment_nt)}
    frag_aln_nogaps_aa = fragments_ch.map{ name, fragment_dir, fragment_nt, fragment_aa -> tuple ("${name}_${fragment_dir.getFileName()}", fragment_aa)}

    frag_raxml_nt(frag_aln_nogaps_nt)
    frag_raxml_aa(frag_aln_nogaps_aa)

    frag_tree_nt = frag_raxml_nt.out
    frag_tree_aa = frag_raxml_aa.out

    // meta map between the original full id and the fragment_id (composed of fasta name and fragment number)
    map_full_frag = fragments_ch.map{ name, fragment_dir, fragment_nt, fragment_aa -> tuple (name, "${name}_${fragment_dir.getFileName()}")}
    map_frag_full = fragments_ch.map{ name, fragment_dir, fragment_nt, fragment_aa -> tuple ("${name}_${fragment_dir.getFileName()}", name)}
    //map_full_frag_small = fragments_ch.map{ name, fragment_dir, fragment_nt, fragment_aa -> tuple (name, "${fragment_dir.getFileName()}")}
    map_all = fragments_ch.map{ name, fragment_dir, fragment_nt, fragment_aa -> tuple (name, "${name}_${fragment_dir.getFileName()}","${fragment_dir.getFileName()}")}
    fragment_names_c = fragments_ch.map{ name, fragment_dir, fragment_nt, fragment_aa -> tuple (name,"${fragment_dir.getFileName()}")}.groupTuple()
    // meta map between the fragment_id (composed of fasta name and fragment number) and the correct original full nt aln file w/o gaps
    frag_corresponding_full_nt_aln_nogaps = map_full_frag.combine(remove_gaps.out.fna, by: 0).map {name, fragment_name, aln -> tuple (fragment_name, aln)}
    // meta map between the fragment_id (composed of fasta name and fragment number) and the correct original full aa aln file w/o gaps
    frag_corresponding_full_aa_aln_nogaps = map_full_frag.combine(remove_gaps.out.faa, by: 0).map {name, fragment_name, aln -> tuple (fragment_name, aln)}
    // meta map between the fragment_id (composed of fasta name and fragment number) and the correct original full raw aa aln file
    frag_corresponding_full_aa_aln_raw = map_full_frag.combine(raw_aln_aa, by: 0).map {name, fragment_name, aln -> tuple (fragment_name, aln)}
    // meta map between the fragment_id (composed of fasta name and fragment number) and the correct original full raw nt aln file
    frag_corresponding_full_nt_aln_raw = map_full_frag.combine(raw_aln_nt, by: 0).map {name, fragment_name, aln -> tuple (fragment_name, aln)}
    // meta map between the fragment_id (composed of fasta name and fragment number) and the correct internal2input file
    frag_internal2input_c = map_full_frag.combine(internal2input_c, by: 0).map {name, fragment_name, tsv -> tuple (fragment_name, tsv)}
    frag_input2interal_c = map_full_frag.combine(input2internal_c, by: 0).map {name, fragment_name, tsv -> tuple (fragment_name, tsv)}


    frag_newicks_ch = 
        frag_raxml2drawing(
            frag_raxml_nt.out.combine(frag_raxml_aa.out, by: 0)
        ).combine(
            frag_corresponding_full_nt_aln_nogaps, by: 0 // simply needed to count the number of taxa for plotting height estimation
        )

    // if outgroups are provided, try reroot the fragment trees as well
    if (params.outgroup != 'NA' ) {
        frag_reroot(frag_newicks_ch.map{name, nt_tree, aa_tree, nt_nogaps_fasta -> tuple (name, nt_tree, aa_tree)})
        frag_tree_nt = frag_reroot.out.nt
        frag_tree_aa = frag_reroot.out.aa
        // update the newicks_ch with the re-rooted trees
        frag_newicks_ch = frag_tree_nt.join(frag_tree_aa, by: 0).combine(frag_corresponding_full_nt_aln_nogaps, by: 0)
    }

    frag_nw_display(frag_newicks_ch)

    // check according to the initial breakpoint array if they are significant
    // ...

    // now run CODEML for each fragment separately
    // fragment nt aln w/o gaps and the corresponding nt tree w/o bootstraps
    frag_aln_tree_ch = frag_aln_nogaps_nt.join(frag_barefoot(frag_raxml_nt.out))

    frag_codeml_combine(
            frag_codeml_run(
                frag_codeml_built(frag_aln_tree_ch).ctl_files.transpose().combine(frag_aln_tree_ch, by: 0)
            ).mlc_files.groupTuple(by: 1, size: 6).map { it -> tuple ( it[0][1], it[1], it[2])}//[bats_mx1_fragment_1, bats_mx1_fragment_1_F1X4, [codeml_F1X4_M1a.mlc, codeml_F1X4_M0.mlc, codeml_F1X4_M8a.mlc, codeml_F1X4_M7.mlc, codeml_F1X4_M8.mlc, codeml_F1X4_M2a.mlc]]
    )//[bats_mx1_fragment_1, F61, /home/martin/git/poseidon/work/7b/43745a2f2434cd51c0ea91945261e6/codeml_F61.all.mlc]    

    breakpoints_ch = map_full_frag.combine(build_fragments.out.breakpoints, by: 0).map {name, fragment_name, bp_tsv -> tuple (fragment_name, bp_tsv)}

    // build LaTeX summary table for each fragment separately
    frag_pdflatex_single(
            frag_tex_built(
                frag_codeml_combine.out
                    .combine(frag_corresponding_full_aa_aln_nogaps, by: 0)
                    .combine(frag_corresponding_full_aa_aln_raw, by: 0)
                    .combine(frag_internal2input_c, by: 0)
                    .combine(breakpoints_ch, by: 0)
            ).tex_files
    )

    /*******************************
    Build combined TeX and PDF, also add fragments */
    fragment_tex_dirs = map_frag_full.combine(
        frag_tex_built.out.tex_dir, by: 0
    ).map{name_frag, name, tex_dir -> tuple (name, tex_dir)}
    .groupTuple(by: 0)

    pdflatex_full(
            tex_combine(
                tex_built.out.tex_dir.groupTuple(by: 0)//.map { name, tex_dirs -> tuple ( name, tex_dirs[0], tex_dirs[1], tex_dirs[2] ) }
                .combine(fragment_tex_dirs, by: 0)
                .join(build_fragments.out.breakpoints, by: 0)
            )
    )//[bats_mx1, full, [bats_mx1_codeml.pdf, bats_mx1_gaps_codeml.pdf]]


    /*******************************
    HTML */

    // build also HTML output for the fragments, if any
    frag_aln_html(
        map_all
        .join(tex_built.out.gap_start2gap_length, by: 0)
        .combine(raw_aln_html, by: 0).map {name, name_fragment, fragment, gap_start2gap_length, html_aln -> tuple (name_fragment, name, fragment, gap_start2gap_length, html_aln)}
        .combine(breakpoints_ch, by: 0)
        .combine(frag_corresponding_full_nt_aln_raw, by: 0)
        .combine(frag_corresponding_full_aa_aln_raw, by: 0)
    )//.all.view()

    frag_aa_bp_with_gaps = frag_aln_html.out.aa_bp_with_gaps
        //.concat(frag_aln_html.out.dummy)
        .unique()
        //.view()
    frag_aa_bp_with_gaps_for_full_aln = map_frag_full.combine(frag_aa_bp_with_gaps, by: 0).map{name_frag, name, csv -> tuple (name, csv)}.unique()

    gap_adjusted_start_end = frag_aln_html.out.gap_adjusted_start_end
    frag_mlc_files_c = frag_codeml_run.out.mlc_files.groupTuple(by: 1, size: 6).map { it -> tuple ( it[0][1], it[2] )}.groupTuple(by: 0).map { it -> tuple (it[0], it[1][0], it[1][1], it[1][2]) }
    html_frag('fragment',
        frag_aln_html.out.all
            .map{name_frag, name, frag, html, aa_aln, nt_aln -> tuple (name_frag, html, aa_aln)}
            //.join(map_full_frag.combine(fasta_input_ch, by: 0).map{name, name_frag, fasta -> tuple (name)})//need this for correct out folder
            .join(frag_codeml_built.out.ctl_dir)
            .join(frag_tree_nt)
            .join(frag_tree_aa)
            .join(map_full_frag.combine(fasta_input_ch, by: 0).map{name, name_frag, fasta -> tuple (name_frag, fasta)})
            .join(frag_internal2input_c)
            .join(frag_input2interal_c)
            .join(map_full_frag.combine(gard_process.out.html, by: 0).map{name, name_frag, gard_html -> tuple (name_frag, gard_html)})
            .join(map_full_frag.combine(model_selection.out.model, by: 0).map{name, name_frag, model -> tuple (name_frag, model)})
            .join(tex_combine.out.map {it -> tuple (it[0], it[2])}//[bats_mx1, [bats_mx1_codeml.tex, bats_mx1_gaps_codeml.tex]]
                .combine(map_full_frag, by: 0).map{name, codeml_tex, name_frag -> tuple(name_frag, codeml_tex)})
            .join(frag_tex_built.out.tex_dir.groupTuple(by: 0))
//            .join(tex_built.out.tex_dir.groupTuple(by: 0)
//                .combine(map_full_frag, by: 0).map{name, tex, name_frag -> tuple(name_frag, tex)})
            .join(frag_mlc_files_c)//[bats_mx1_fragment_1, [codeml_F1X4_M0.mlc, codeml_F1X4_M8a.mlc, codeml_F1X4_M1a.mlc, codeml_F1X4_M8.mlc, codeml_F1X4_M2a.mlc, codeml_F1X4_M7.mlc], [codeml_F3X4_M0.mlc, codeml_F3X4_M1a.mlc, codeml_F3X4_M7.mlc, codeml_F3X4_M2a.mlc, codeml_F3X4_M8.mlc, codeml_F3X4_M8a.mlc], [codeml_F61_M7.mlc, codeml_F61_M8.mlc, codeml_F61_M1a.mlc, codeml_F61_M2a.mlc, codeml_F61_M0.mlc, codeml_F61_M8a.mlc]]
            .join(frag_tex_built.out.lrt_params.groupTuple(by: 0))//[bats_mx1_fragment_1, [F3X4.lrt, F1X4.lrt, F61.lrt]]
//            .join(tex_built.out.lrt_params.groupTuple(by: 0)//[bats_mx1_fragment_1, [F3X4.lrt, F1X4.lrt, F61.lrt]]
//                .combine(map_full_frag, by: 0).map{name, lrt, name_frag -> tuple(name_frag, lrt)})
            .join(pdflatex_full.out.map {it -> tuple (it[0], it[2])}
                .combine(map_full_frag, by: 0).map{name, pdf, name_frag -> tuple(name_frag, pdf)})
            .join(frag_aln_html.out.all.map{name_frag, name, frag, html, aa_aln, nt_aln -> tuple (name_frag, nt_aln)})
            .join(frag_nw_display.out.map { it -> tuple (it[0], it[1], it[2], it[3]) })
            .join(map_full_frag.combine(model_selection.out.log, by: 0).map{name, name_frag, log -> tuple (name_frag, log)})
            .join(translatorx.out.aa
                .combine(map_full_frag, by: 0).map{name, aa, name_frag -> tuple(name_frag, aa)})
            .join(remove_gaps.out.fna
                .combine(map_full_frag, by: 0).map{name, fna, name_frag -> tuple(name_frag, fna)})
            .join(remove_gaps.out.faa
                .combine(map_full_frag, by: 0).map{name, faa, name_frag -> tuple(name_frag, faa)})
            .join(frag_aa_bp_with_gaps, by: 0)
            .join(gap_adjusted_start_end, by: 0)//a csv with the gap-adjusted start and end pos for a certain fragment
            .join(gard_process.out.recombination
                .combine(map_full_frag, by: 0).map{name, recombination, name_frag -> tuple(name_frag, recombination)})//[bats_mx1_fragment_1, true]
            //.view()
    )

    // Frag CODEML SUMMARY
    frag_html_codeml('fragment',
                html_frag.out.index
                .join(html_frag.out.tex_dir)
                .join(frag_tex_built.out.lrt_params.groupTuple(by: 0))//[bats_mx1_fragment_1, [F3X4.lrt, F1X4.lrt, F61.lrt]]
                .join(map_full_frag.combine(fragment_names_c, by: 0).map{name, frag_name, frag -> tuple (frag_name, frag)}, by: 0)
                //.view()
    )
    ///////////////////////////////////////////////////////////////////////////////////// Fragments ending

    // MAIN SUMMARY

    mlc_files_c = codeml_run.out.mlc_files.groupTuple(by: 1, size: 6).map { it -> tuple ( it[0][1], it[2] )}.groupTuple(by: 0).map { it -> tuple (it[0], it[1][0], it[1][1], it[1][2]) }
    html('full_aln', 
        checked_aln_html
            .join(checked_aln_aa)
            //.join(checked_aln_aa.map{name, faa -> name})//need this for correct out folder, bc/ fragments above need this
            .join(codeml_built.out.ctl_dir)
            .join(tree_nt)
            .join(tree_aa)
            .join(fasta_input_ch)
            .join(internal2input_c)
            .join(input2internal_c)
            .join(gard_process.out.html)
            .join(model_selection.out.model)
            .join(tex_combine.out.map {it -> tuple (it[0], it[2]) })//[bats_mx1, [bats_mx1_codeml.tex, bats_mx1_gaps_codeml.tex]]
            .join(tex_built.out.tex_dir.groupTuple(by: 0))
            .join(mlc_files_c)//[bats_mx1, [codeml_F1X4_M0.mlc, codeml_F1X4_M8a.mlc, codeml_F1X4_M1a.mlc, codeml_F1X4_M8.mlc, codeml_F1X4_M2a.mlc, codeml_F1X4_M7.mlc], [codeml_F3X4_M0.mlc, codeml_F3X4_M1a.mlc, codeml_F3X4_M7.mlc, codeml_F3X4_M2a.mlc, codeml_F3X4_M8.mlc, codeml_F3X4_M8a.mlc], [codeml_F61_M7.mlc, codeml_F61_M8.mlc, codeml_F61_M1a.mlc, codeml_F61_M2a.mlc, codeml_F61_M0.mlc, codeml_F61_M8a.mlc]]
            .join(tex_built.out.lrt_params.groupTuple(by: 0))//[bats_mx1, [F3X4.lrt, F1X4.lrt, F61.lrt]]
            .join(pdflatex_full.out.map {it -> tuple (it[0], it[2])})
            .join(checked_aln_nt)
            .join(nw_display.out.map { it -> tuple (it[0], it[1], it[2], it[3]) })
            .join(model_selection.out.log)
            .join(translatorx.out.aa)
            .join(remove_gaps.out.fna)
            .join(remove_gaps.out.faa)
            .join(frag_aa_bp_with_gaps_for_full_aln)//aln_length_w_gaps default for full aln
            .join(gap_adjusted_start_end_main)
            .join(gard_process.out.recombination)//[bats_mx1, true]
            //.view()
    )

    // CODEML SUMMARY
    html_codeml('full_aln',
                html.out.index
                .join(html.out.tex_dir)
                .join(tex_built.out.lrt_params.groupTuple(by: 0))//[bats_mx1, [F3X4.lrt, F1X4.lrt, F61.lrt]]
                .join(fragment_names_c, by: 0)
    )

    // PARAMETER SUMMARY
    html_params('full_aln', html.out.index)

    // if fragments!!!
    // Recombination summary HTML
    html_recomb(
        html.out.index
        .join(gard_process.out.html)
        .join(tree_nt)
        .join(map_frag_full.combine(frag_tree_nt, by: 0).map{frag_name, name, trees -> tuple (name, trees)}.groupTuple(by:0))
        //.view()
    )    

    // Summarize all output files


}


/*************  
* --help
*************/

def helpMSG() {
    c_green = "\033[0;32m";
    c_reset = "\033[0m";
    c_yellow = "\033[0;33m";
    c_blue = "\033[0;34m";
    c_dim = "\033[2m";
    log.info """
    ____________________________________________________________________________________________
    
    PoSeiDon -- Positive Selection Detection and Recombination Analysis

    ${c_yellow}Usage example:${c_reset}
    nextflow run poseidon.nf --fasta '*/*.fasta' 

    ${c_yellow}Input:${c_reset}
    ${c_green} --fasta ${c_reset}       '*.fasta'           -> one FASTA file per transcriptome assembly
    ${c_dim}  ..change above input to csv:${c_reset} ${c_green}--list ${c_reset}
    
    ${c_yellow}General options:${c_reset}
    --cores             max cores for local use [default: $params.cores]
    --memory            memory limitations for polisher tools in GB [default: $params.memory]
    --output            name of the result folder [default: $params.output]
    --reference         resulting amino acid changes and sites will be reported according to this species (FASTA id) [default: $params.reference]
    --root              outgroup species (FASTA id) for tree rooting; comma-separated [default: $params.outgroup]
    --bootstrap         number of bootstrap calculations [default: $params.bootstrap]

    ${c_yellow}Model parameters:${c_reset}
    --model             nucleotide model used for recombination analysis [default: $params.model]
    --model_rc          model rate classes [default: $params.model_rate_classes]
    --model_sm          model selection method [default: $params.model_selection_method] 
    --model_rl          model rejection level [default: $params.model_rejection_level]

    ${c_yellow}Recombination parameters (GARD):${c_reset}
    --gard_rv           GARD rate variation [default: $params.gard_rate_variation]
    --gard_rc           GARD rate classes [default: $params.gard_rate_classes]
    --kh                use insignificant breakpoints (based on KH test) for fragment calcuations [default: $params.kh]

    ${c_dim}Nextflow options:
    -with-report rep.html    cpu / ram usage (may cause errors)
    -with-dag chart.html     generates a flowchart for the process tree
    -with-timeline time.html timeline (may cause errors)

    ${c_yellow}LSF computing:${c_reset}
    For execution of the workflow on a HPC with LSF adjust the following parameters:
    --databases         defines the path where databases are stored [default: $params.cloudDatabase]
    --workdir           defines the path where nextflow writes tmp files [default: $params.workdir]
    --cachedir          defines the path where images (singularity) are cached [default: $params.cachedir] 

    Profile:
    Merge profiles comma-separated
    -profile                 local,docker
                             local,conda
                             lsf,docker,singularity (adjust workdir and cachedir according to your HPC config)
                             slurm,conda (adjust workdir and cachedir according to your HPC config)
                             gcloud,docker (GCP google-lifescience with docker)
                             ${c_reset}
    """.stripIndent()
}

