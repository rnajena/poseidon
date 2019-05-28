# PoSeiDon

Here we present __PoSeiDon__, a pipeline to detect significant positively selected sites and possible recombination events in analignment of multiple coding sequences. Sites that undergo positive selection can give you insights in the evolutionary history of your sequences, for example showing you important mutation hot spots, accumulated as results of virus-host arms races during evolution.

We provide all ruby scripts needed to run the PoSeiDon pipeline.

__Please note__: we aimed that with these scripts the pipeline can be run _out
of the box_, however, PoSeiDon relies on a variety of different third-party
tools (see below). Binaries for most tools are also included in this repository
(`tools`) and PoSeiDon assumes them to be located in this folder. The larger
software package for HYPHY can be downloaded here directly and needs to be added
and extracted manually to the `tools` folder:

* <a href="https://www.rna.uni-jena.de/supplements/poseidon/hyphy.zip">hyphy.zip</a>
<!--* <a href="https://www.rna.uni-jena.de/supplements/poseidon/openmpi.zip">openmpi.zip</a>-->

Furthermore, you will need inkscape, pdflatex, ruby (tested with v2.4.2) and
some ruby gems (packages) as well as mpirun (Open MPI; tested with v2.0.2). If
you don't have anything of this installed, you can try on a Linux system:

````
apt-get install ruby
gem install bio
gem install mail
gem install encrypted_strings

apt-get install inkscape

apt-get install texlive-latex-base

apt-get install openmpi-bin

apt-get install hyphy-mpi
````

__We heavily recommend__ to use our Docker image that can be easily executed without the need to install tools manually.:

````
docker run mhoelzer/poseidon <TODO>
````

## Workflow of the PoSeiDon pipeline and example output

<a target="_blank" href="https://github.com/hoelzer/poseidon/blob/master/images/pipeline_landscape.pdf"><img src="https://github.com/hoelzer/poseidon/blob/master/images/pipeline_landscape.png" alt="PoSeiDon workflow" /></a>

The PoSeiDon pipeline comprises in-frame alignment of homologous protein-coding sequences, detection of putative recombination events and evolutionary breakpoints, phylogenetic reconstructions and detection of positively selected sites in the full alignment and all possible fragments. Finally, all results are combined and visualized in a user-friendly and clear HTML web page. The resulting alignment fragments are indicated with colored bars in the HTML output.

Please find an example output of the pipeline <a href="http://www.rna.uni-jena.de/supplements/mx1_bats/full_aln/">here</a>. (<a href="https://doi.org/10.1128/JVI.00361-17">Fuchs _et al_., 2017, Journal of Virology</a>)

### The PoSeiDon pipeline is based on the following tools and scripts:

* TranslatorX (v1.1), Abascal et al. (2010);  <a target="_blank" href="https://www.ncbi.nlm.nih.gov/pubmed/20435676">20435676</a>
* Muscle (v3.8.31), Edgar (2004);  <a target="_blank" href="https://www.ncbi.nlm.nih.gov/pubmed/15034147">15034147</a>
* RAxML (v8.0.25), Stamatakis (2014);  <a target="_blank" href="https://www.ncbi.nlm.nih.gov/pubmed/24451623">24451623</a>
* Newick Utilities (v1.6), Junier and Zdobnov (2010);  <a target="_blank" href="https://www.ncbi.nlm.nih.gov/pubmed/20472542">20472542</a>
* MODELTEST , Posada and Crandall (1998);  <a target="_blank" href="https://www.ncbi.nlm.nih.gov/pubmed/9918953">9918953</a>
* HyPhy (v2.2), Pond et al. (2005);  <a target="_blank" href="https://www.ncbi.nlm.nih.gov/pubmed/15509596">15509596</a>
* GARD , Pond et al. (2006);  <a target="_blank" href="https://www.ncbi.nlm.nih.gov/pubmed/17110367">17110367</a>
* PaML/CodeML (v4.8), Yang (2007);  <a target="_blank" href="https://www.ncbi.nlm.nih.gov/pubmed/17483113">17483113</a>
* Ruby (v2.3.1)
* Inkscape (v0.48.5)
* pdfTeX (v3.14) 

## Parameters

Most of the PoSeiDon parameters are optional and are explained here in detail.

### Input
Mandatory. Your input FASTA file must follow the format:

````
>Myotis_lucifugus Mx1 Gene
ATGGCGATCGAGATACGATACGTA...
>Myotis_davidii Mx1 Gene
ATGGCGGTCGAGATAAGATACGTT...
````

All sequences must have a correct open reading frame, are only allowed to contain nucleotide characters [A|C|G|T] and no internal stop codon.

Sequence IDs must be unique until the first occurrence of a space.

### Reference
Optional. Default: use first sequence ID as reference. You can define <b>one</b> species ID from your multiple FASTA file as a reference species. Positively selected sites and corresponding amino acids will be drawn in respect to this species. The ID must match the FASTA header until the occurence of the first space. For example, if you want <i>Myotis lucifugus</i> as your reference species and your FASTA file contains:

````>Myotis_lucifugus Mx1 Gene
ATGGCGATCGAGATACGATACGTA...
````

use

````Myotis_lucifugus````

as parameter to set the reference species. Per default the first ID occurring in the multiple FASTA file will be used.

### Outgroup
Optional. Default: trees are unrooted. You can define <b>one</b> or <b>multiple</b> (comma separated) species IDs as outgroup. All phylogenetic trees will be rooted according to this species. For example, if your multiple FASTA file contains

````>Myotis_lucifugus Mx1 Gene
ATGGCGATCGAGATACGATACGTA...
>Myotis_davidii Mx1 Gene
ATGGCGGTCGAGATAAGATACGTT...
>Pteropus_vampyrus Mx1 Gene
ATGGCCGTAGAGATTAGATACTTT...
>Eidolon_helvum Mx1 Gene
ATGCCCGTAGAGAATAGATACTTT...
````

you can define:

````Pteropus_vampyrus,Eidolon_helvum````

to root all trees in relation to this two species.

### Use also insignificant breakpoints
Optional. Default: false. With this parameter you can decide if insignificant breakpoints should be taken into account. All breakpoints are tested for significant topological incongruence using a Kashino Hasegawa (KH) test [Kishino, H. and Hasegawa, M. (1989)]. KH-insignificant breakpoints most frequently arise from variation in branch lengths between segments. Nevertheless, taking KH-insignificant breakpoints into account could be interesting, because we already observed putative positively selected sites in fragments without any significant topological incongruence. KH-insignificant fragments are marked in the final output, as they might not occur from real recombination events.

Per default only significant breakpoints are used for further calculations.

Please also keep in mind that using also insignificant breakpoints can extend the run time of PoSeiDon from minutes to hours, depending on the number of detected breakpoints.

### Use your own parameters

Currently, we don't provide full access to the parameters used within PoSeiDon through the web interface __[the web serice is currently under maintenance due to web page changes]__. In a future release, we will provide a local version of the pipeline for download including full access to the parameter settings of all executed tools. If you want to change parameters (e.g. for RAxML) now, just run the pipeline and PoSeiDon will also generate a 'Parameters' sub page (like <a href="https://www.rna.uni-jena.de/supplements/mx1_bats/full_aln/params.html">this</a>) in the final output, allowing access to all executed commands. With this, certain parts of the pipeline can be rerun locally using the provided commands and output files.
