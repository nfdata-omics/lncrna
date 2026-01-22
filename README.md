<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nfdata-omics-lncrna_logo.png">
    <img alt="nf-core/lncrna" src="docs/images/nfdata-omics-lncrna_logo.png">
  </picture>
</h1>

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/nfdata-omics/lncrna)
[![GitHub Actions CI Status](https://github.com/nfdata-omics/lncrna/actions/workflows/nf-test.yml/badge.svg)](https://github.com/nfdata-omics/lncrna/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/nfdata-omics/lncrna/actions/workflows/linting.yml/badge.svg)](https://github.com/nfdata-omics/lncrna/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.5.1-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.5.1)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nfdata-omics/lncrna)

## Introduction

**nfdata-omics/lncrna** is a bioinformatics pipeline that processes RNA‑seq data to quantify known lncRNAs and optionally discover novel lncRNAs. It performs FASTQ QC and alignment, transcript assembly and novelty filtering, coding‑potential assessment, genomic‑context classification, and reporting. By default it runs in known‑only-lncrnas mode; enabling the discovery branch adds identification, coding-potential assessment and a final evaluation of lncRNAs.

![nfdata-omics/lncrna metro map](docs/images/nfdata-omics-lncrna-pipeline.png)

- FASTQ QC and trimming; optional rRNA removal and filtering
- Alignment and quantification (STAR or HISAT2)
- Transcript assembly (StringTie) and merge
- Pseudo-allignment with Salmon.
- Optional novel lncRNAs identification: Coding-potential (CPAT/FEELnc/PLEK), filter against known lncRNAs, classification.
- Final HTML report and aggregated QC (MultiQC)

<!-- Add a tube map or workflow figure here if desired -->

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set up Nextflow. Test your setup with `-profile test` before running on actual data.

Prepare a samplesheet CSV:

```csv
sample,fastq_1,fastq_2
CTRL_1,/data/CTRL_1_R1.fastq.gz,/data/CTRL_1_R2.fastq.gz
```

Run in known‑only mode (default):

```bash
nextflow run main.nf \
  --input assets/samplesheet.csv \
  --fasta /path/to/genome.fa \
  --gtf /path/to/genes.gtf \
  --outdir /path/to/output \
  -profile conda \
  -resume
```

Enable novel lncRNA discovery:

```bash
nextflow run main.nf \
  --input assets/samplesheet.csv \
  --fasta /path/to/genome.fa \
  --gtf /path/to/genes.gtf \
  --outdir /path/to/output_novel \
  --novel_lncrnas true \
  -profile conda \
  -resume
```

> [!WARNING]
> Provide pipeline parameters via the CLI or Nextflow `-params-file`. Custom config files (`-c`) can adjust executor and resource configuration but should not define parameters; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Outputs

### Known lncRNA assessment

- lncrna_gtf_filtering/by_classcode: filtered GTF and classcode stats
- lncrna_final_annotation: final combined annotation and splits
- expression/quantification: per‑sample quantification outputs
- multiqc: aggregated QC report

### Novel lncRNAs discovery

- lncrna_transcript_filtering/transcripts_length: length‑filtered outputs
- lncrna_transcript_filtering/transcripts_exons: exon‑filtered outputs
- lncrna_prediction/cpat: CPAT predictions
- lncrna_prediction/cpat/models: CPAT models when built
- lncrna_prediction/plek and lncrna_prediction/feelnc: tool‑specific outputs
- lncrna_prediction/combined_predictions: consensus results
- expression/matrix: merged count matrix and gene metadata
- lncrna_report: final HTML report and summary text
- multiqc: aggregated QC report

## Credits

nfdata-omics/lncrna was originally written by K. Ruiz-Ceja, M. Bonfanti, ....

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nfdata-omics/lncrna for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
