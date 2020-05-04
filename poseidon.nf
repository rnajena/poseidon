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
include model_selection from './modules/model_selection'
include {codeml_run; codeml_built; codeml_combine} from './modules/codeml'
include {tex_built; pdflatex} from './modules/tex'

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

    // check input FASTA
    check_fasta_format(fasta_input_ch)
    
    // TODO rndm subsampling if too many sequences, PAML is not inteded for >100 sequences

    // align, check alignment, and remove gaps
    remove_gaps(
        check_aln(
            translatorx(
                check_fasta_format.out.fasta).join(check_fasta_format.out.log)))

    // raw alignments from translatorx
    raw_aln_aa = translatorx.out.map {name, nt_aln, aa_aln, html -> [name, aa_aln] }
    raw_aln_nt = translatorx.out.map {name, nt_aln, aa_aln, html -> [name, nt_aln] }

    // ID mappings
    internal2input_c = check_fasta_format.out.internal2input
    input2internal_c = check_fasta_format.out.input2internal

    // build nt and aa phylogeny
    raxml_nt(remove_gaps.out.fna)
    raxml_aa(remove_gaps.out.faa)

    // convert RAxML output newicks for drawing
    newicks_ch = 
        raxml2drawing(
            raxml_nt.out.concat(raxml_aa.out)
        ).combine(
            remove_gaps.out.fna, by: 0 // simply needed to count the number of taxa for plotting height estimation
        )

    // if outgroups are provided, try reroot the trees
    // TODO, not implemented yet
    if (params.root != 'NA' ) {
        reroot(newicks_ch)
    }

    // draw trees
    nw_display(newicks_ch)

    // full nt aln w/o gaps and the corresponding nt tree w/o bootstraps
    aln_tree_ch = remove_gaps.out.fna.join(barefoot(raxml_nt.out))

    // model selection
    model_selection(aln_tree_ch)

    // GARD breakpoint analysis -- TODO: I could skip this for now and just do everything on the full aln
    //gard(remove_gaps.out.fna.join(model_selection.out))
    
    // CODEML
    //codeml_combine(
            codeml_run(
                codeml_built(aln_tree_ch).transpose().combine(aln_tree_ch, by: 0)
            ).view()//.groupTuple(by: 1, size: 6).map { it -> tuple ( it[0][1], it[1], it[2])} //[bats_mx1, F1X4, [/home/martin/git/poseidon/work/a8/0b99e44e367ab9c65191ace42a0370/codeml_F1X4_M1a.mlc, /home/martin/git/poseidon/work/c9/a7cd64899060251d3242c6e0860a2e/codeml_F1X4_M0.mlc, /home/martin/git/poseidon/work/b2/06b79f4f07a978db0efd7aa5db9c86/codeml_F1X4_M8a.mlc, /home/martin/git/poseidon/work/18/52f0dabb8c9793de55fabf57df2558/codeml_F1X4_M7.mlc, /home/martin/git/poseidon/work/a8/752ef1b9a3545d1068f0afd8cd4d36/codeml_F1X4_M8.mlc, /home/martin/git/poseidon/work/d4/f5710b4dd4f9517861a2f8db522df5/codeml_F1X4_M2a.mlc]]
    //) //[bats_mx1, F61, /home/martin/git/poseidon/work/7b/43745a2f2434cd51c0ea91945261e6/codeml_F61.all.mlc]


    //codeml_combine.out.view()
    // LaTeX summary table
    //pdflatex(
    //    tex_built(
    //        codeml_combine.out.combine(remove_gaps.out.faa, by: 0).combine(raw_aln_aa, by: 0).combine(internal2input_c, by: 0)
    //    )
    //)


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
    
    ${c_yellow}Options:${c_reset}
    --cores             max cores for local use [default: $params.cores]
    --memory            memory limitations for polisher tools in GB [default: $params.memory]
    --output            name of the result folder [default: $params.output]
    --reference         resulting positions will be reported according to this species [default: $params.reference]
    --root              outgroup species for tree rooting; comma-separated [default: $params.root]
    --bootstrap         number of bootstrap calculations [default: $params.bootstrap]

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

