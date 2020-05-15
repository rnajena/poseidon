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

include {raxml_nt; raxml_aa; raxml2drawing; nw_display; barefoot} from './modules/tree'
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

include {html_main; html_codeml; html_params; frag_aln_html} from './modules/html'

include build_fragments from './modules/fragment'

/************************** 
* DATABASES
**************************/

/* Comment section:
The Database Section is designed to "auto-get" pre prepared databases.
It is written for local use and cloud use via params.cloudProcess.
*/


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
    raxml_nt(remove_gaps.out.fna)
    raxml_aa(remove_gaps.out.faa)

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
    // TODO, not implemented yet
    if (params.root != 'NA' ) {
        reroot(newicks_ch)
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

    /*******************************
    -- FRAGMENTS --
    If breakpoints ... */ 
    fragments_ch = build_fragments(remove_gaps.out.fna.join(gard_process.out.bp).join(gard_process.out.recombination)).fragments.transpose()

    // for each fragment, build a new unrooted tree for CODEML
    frag_aln_nogaps_nt = fragments_ch.map{ name, fragment_dir, fragment_nt, fragment_aa -> tuple ("${name}_${fragment_dir.getFileName()}", fragment_nt)}
    frag_aln_nogaps_aa = fragments_ch.map{ name, fragment_dir, fragment_nt, fragment_aa -> tuple ("${name}_${fragment_dir.getFileName()}", fragment_aa)}

    frag_raxml_nt(frag_aln_nogaps_nt)
    frag_raxml_aa(frag_aln_nogaps_aa)

    // meta map between the original full id and the fragment_id (composed of fasta name and fragment number)
    map_full_frag = fragments_ch.map{ name, fragment_dir, fragment_nt, fragment_aa -> tuple (name, "${name}_${fragment_dir.getFileName()}")}
    map_full_frag_small = fragments_ch.map{ name, fragment_dir, fragment_nt, fragment_aa -> tuple (name, "${fragment_dir.getFileName()}")}
    map_all = fragments_ch.map{ name, fragment_dir, fragment_nt, fragment_aa -> tuple (name, "${name}_${fragment_dir.getFileName()}","${fragment_dir.getFileName()}")}
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

    frag_newicks_ch = 
        frag_raxml2drawing(
            frag_raxml_nt.out.combine(frag_raxml_aa.out, by: 0)
        ).combine(
            frag_corresponding_full_nt_aln_nogaps, by: 0 // simply needed to count the number of taxa for plotting height estimation
        )
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
    Build combined TeX and PDF */
    pdflatex_full(
            tex_combine(tex_built.out.tex_files.groupTuple(by: 0).map { it -> tuple ( it[0], it[2][0], it[2][1], it[2][2] ) })
    )//[bats_mx1, full, [/home/martin/git/poseidon/work/36/53218cb7db8f91aa7243a6216249f1/bats_mx1_codeml.pdf, /home/martin/git/poseidon/work/36/53218cb7db8f91aa7243a6216249f1/bats_mx1_gaps_codeml.pdf]]

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
    ).view()

    // MAIN SUMMARY
    // get correct trees
    if (params.root != 'NA') {
        tree_nt = reroot.out.nt
        tree_aa = reroot.out.aa
    } else {
        tree_nt = raxml_nt.out
        tree_aa = raxml_aa.out
    }

    mlc_files_c = codeml_run.out.mlc_files.groupTuple(by: 1, size: 6).map { it -> tuple ( it[0][1], it[2] )}.groupTuple(by: 0).map { it -> tuple (it[0], it[1][0], it[1][1], it[1][2]) }
    html_main('full_aln', 
        checked_aln_html
            .join(checked_aln_aa)
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
            //.view()
    )




    // CODEML SUMMARY
    //fragment_names_c = Channel.empty()
    html_codeml('full_aln',
                html_main.out.index
                .join(html_main.out.tex_dir)
                .join(tex_built.out.lrt_params.groupTuple(by: 0))//[bats_mx1, [F3X4.lrt, F1X4.lrt, F61.lrt]]
                //fragment_names_c
    )

    // PARAMETER SUMMARY
    html_params('full_aln', html_main.out.index)

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
    --reference         resulting positions will be reported according to this species [default: $params.reference]
    --root              outgroup species for tree rooting; comma-separated [default: $params.root]
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
                             gcloudMartin,docker (GCP google-lifescience with docker)
                             ${c_reset}
    """.stripIndent()
}

