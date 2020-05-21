# PoSeiDon

<img align="right" width="120px" src="https://github.com/hoelzer/poseidon/blob/master/images/poseidon_logo.png" alt="PoSeiDon logo" /> 

![](https://img.shields.io/badge/nextflow-20.01.0-brightgreen)
![](https://img.shields.io/badge/uses-docker-blue.svg)
![](https://img.shields.io/badge/licence-GPL--3.0-lightgrey.svg)

[![Twitter Follow](https://img.shields.io/twitter/follow/martinhoelzer.svg?style=social)](https://twitter.com/martinhoelzer) 

__Please note that the code of PoSeiDon was just transfered to Nextflow so there might be still some bugs. Please feel free to report issues!__

Here we present [PoSeiDon](https://www.biorxiv.org/content/10.1101/2020.05.18.102731v1), a pipeline to detect significant positively selected sites and possible recombination events in an alignment of multiple protein-coding sequences. Sites that undergo positive selection provide insights in the evolutionary history of your sequences, for example showing important mutation hot spots, accumulated as results of virus-host _arms races_ during evolution.

PoSeiDon relies on a variety of different third-party tools (see below). But don't worry, we encapsulated each tool in its own [Docker](https://www.docker.com/resources/what-container) container and connected them in the Workflow Management System [Nextflow](https://www.nextflow.io/). 

## Installation

You only need Nextflow (version 20.+) and Docker installed to run the pipeline. All dependencies will be pulled automatically. 

Either run PoSeiDon by cloning this repository:
```bash
git clone https://github.com/hoelzer/poseidon.git
cd poseidon
nextflow run poseidon.nf --help
```

or let Nextflow do the pull
```bash
nextflow run hoelzer/poseidon --help
```

We recommend using a specific release of PoSeiDon (as soon as their are releases) via 

```bash
nextflow run hoelzer/poseidon -r v0.1 --help
```

### Update

Depending on you installation procedure, update the pipeline via `git pull` or `nextflow pull hoelzer/poseidon`. 

## Run

__Important:__ PoSeiDon needs nucleotide sequences with a correct open reading frame as input. In addition, the results heavily depend on your selection of sequences, thus, you might consider running the pipeline multiple times with different samplings of your input sequences. Also the pipeline can't work with too many sequences because in its core PoSeiDon uses CODEML from the [PAML](http://abacus.gene.ucl.ac.uk/software/paml.html) suite that is not inteded for >100 sequences. __Please find a detailed description of the input parameters and settings below.__

### Profiles

Nextflow can be easily executed on different environments like your local machine, a high-performance cluster or the cloud. Different `-profile` are used to tell Nextflow which system should be used. For local execution `-profile local,docker` should be used. More profiles will be added.   

### Examples 

Now let's assume you used Nextflow to pull the PoSeiDon code. 

```bash
# show help 
nextflow run hoelzer/poseidon --help

# run small example on a local machine with (first time this will need some more time because the Docker containers are downloaded)
nextflow run hoelzer/poseidon --fasta ~/.nextflow/assets/hoelzer/poseidon/test_data/bats_mx1_small.fasta \
--cores 4 -profile local,docker

# resume a broken run
nextflow run hoelzer/poseidon --fasta ~/.nextflow/assets/hoelzer/poseidon/test_data/bats_mx1_small.fasta \
--cores 4 -profile local,docker -resume
```

To reproduce the [positive selection results](http://www.rna.uni-jena.de/supplements/mx1_bats/full_aln/) reported in [Fuchs _et al_. (2017), Journal of Virology](https://doi.org/10.1128/JVI.00361-17) run:
```bash
nextflow run hoelzer/poseidon --fasta ~/.nextflow/assets/hoelzer/poseidon/test_data/bats_mx1.fasta \
--cores 4 --kh --outgroup "Pteropus_alecto,Eidolon_helvum,Rousettus_aegyptiacus,Hypsignatus_monstrosus" \
--reference "Myotis_daubentonii" -profile local,docker
```

## Workflow of the PoSeiDon pipeline and example output

<a target="_blank" href="https://github.com/hoelzer/poseidon/blob/master/images/pipeline_landscape.pdf"><img src="https://github.com/hoelzer/poseidon/blob/master/images/pipeline_landscape.png" alt="PoSeiDon workflow" /></a>

The PoSeiDon pipeline comprises in-frame alignment of homologous protein-coding sequences, detection of putative recombination events and evolutionary breakpoints, phylogenetic reconstructions and detection of positively selected sites in the full alignment and all possible fragments. Finally, all results are combined and visualized in a user-friendly and clear HTML web page. The resulting alignment fragments are indicated with colored bars in the HTML output.

### The PoSeiDon pipeline is based on the following tools and scripts:

* TranslatorX (v1.1), Abascal _et al_. (2010);  <a target="_blank" href="https://www.ncbi.nlm.nih.gov/pubmed/20435676">20435676</a>
* Muscle (v3.8.31), Edgar (2004);  <a target="_blank" href="https://www.ncbi.nlm.nih.gov/pubmed/15034147">15034147</a>
* RAxML (v8.0.25), Stamatakis (2014);  <a target="_blank" href="https://www.ncbi.nlm.nih.gov/pubmed/24451623">24451623</a>
* Newick Utilities (v1.6), Junier and Zdobnov (2010);  <a target="_blank" href="https://www.ncbi.nlm.nih.gov/pubmed/20472542">20472542</a>
* MODELTEST, Posada and Crandall (1998);  <a target="_blank" href="https://www.ncbi.nlm.nih.gov/pubmed/9918953">9918953</a>
* HyPhy (v2.2), Pond _et al_. (2005);  <a target="_blank" href="https://www.ncbi.nlm.nih.gov/pubmed/15509596">15509596</a>
* GARD, Pond et al. (2006);  <a target="_blank" href="https://www.ncbi.nlm.nih.gov/pubmed/17110367">17110367</a>
* PaML/CodeML (v4.8), Yang (2007);  <a target="_blank" href="https://www.ncbi.nlm.nih.gov/pubmed/17483113">17483113</a>
* Ruby (v2.3.1)
* Inkscape (v1.0)
* pdfTeX (v3.14) 

## Parameters

Most of the PoSeiDon parameters are optional and are explained below in detail.

### Input

* `--fasta`

Mandatory. Your input FASTA file must follow the format:

```
>Myotis_lucifugus Mx1 Gene
ATGGCGATCGAGATACGATACGTA...
>Myotis_davidii Mx1 Gene
ATGGCGGTCGAGATAAGATACGTT...
```

All sequences must have a correct open reading frame, are only allowed to contain nucleotide characters [A|C|G|T] and no internal stop codon.

Sequence IDs must be unique until the first occurrence of a space.

### Reference

* `--reference`

Optional. Default: use first sequence ID as reference. You can define <b>one</b> species ID from your multiple FASTA file as a reference species. Positively selected sites and corresponding amino acids will be drawn in respect to this species. The ID must match the FASTA header until the occurence of the first space. For example, if you want <i>Myotis lucifugus</i> as your reference species and your FASTA file contains:

```>Myotis_lucifugus Mx1 Gene
ATGGCGATCGAGATACGATACGTA...
```

use

`--reference "Myotis_lucifugus"`

as parameter to set the reference species. Per default the first ID occurring in the multiple FASTA file will be used.

### Outgroup(s)

* `--outgroup`

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

```--outgroup "Pteropus_vampyrus,Eidolon_helvum"```

to root all trees in relation to this two species.

### Use also insignificant breakpoints

* `--kh`

Optional. Default: false. With this parameter you can decide if insignificant breakpoints should be taken into account. All breakpoints are tested for significant topological incongruence using a Kashino Hasegawa (KH) test [Kishino, H. and Hasegawa, M. (1989)](https://link.springer.com/article/10.1007/BF02100115). KH-insignificant breakpoints most frequently arise from variation in branch lengths between segments. Nevertheless, taking KH-insignificant breakpoints into account could be interesting, because we already observed putative positively selected sites in fragments without any significant topological incongruence. KH-insignificant fragments are marked in the final output, as they might not occur from real recombination events.

Per default only significant breakpoints are used for further calculations.

Please also keep in mind that using also insignificant breakpoints can extend the run time of PoSeiDon from minutes to hours, depending on the number of detected breakpoints.

### Other parameters

Please see the `--help` for other parameters (GARD, RAxML, ...) and let us know if you need more customization!

## Citation

If PoSeiDon helps you please cite:

[Martin Hölzer and Manja Marz, "PoSeiDon: a Nextflow pipeline for the detection of evolutionary recombination events and positive selection", _BioRxiv_ (2020)](https://www.biorxiv.org/content/10.1101/2020.05.18.102731v1)
