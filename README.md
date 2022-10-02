# ![bigbio/pride-molecules-indexer](docs/images/Nextflow_logo.png)

**PRIDE Archive nextflow molecules indexer pipeline for complete submission submissions**.

[![install with conda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/)

## Introduction

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

## Quick Start

1. Install [`nextflow`](https://nf-co.re/usage/installation) (`>=20.04.0`)

2. Te tool works with Conda 

3. Start running your own analysis!

    ```bash
     nextflow run main.nf -c nextflow.config -profile ebitools --project_accession PXD004016
    ```

See [usage docs](https://nf-co.re/proteomicslfq/usage) for all the available options when running the pipeline. Or configure the pipeline via
[nf-core launch](https://nf-co.re/launch) from the web or the command line.

## Pipeline Summary

By default, the pipeline currently performs the following:

* Get the result files from a Complete Submission
* Get the releated (spectra) files for each result files 
* Generate Json files for Spectra, PSM and Proteins 

## Documentation

## Credits

bigbio/pride-molecules-indexer was originally written by Yasset Perez-Riverol.

## Contributions and Support

[Discussions & Chat & QA](https://github.com/bigbio/pride-molecule-indexer/discussions)
