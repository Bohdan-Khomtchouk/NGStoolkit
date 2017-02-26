# NGStoolkit

## About

`NGStoolkit` is a one-stop shop next-generation sequencing analysis toolkit in Bash that allows you to run a complete NGS analysis pipeline (e.g., RNA-seq) by issuing just one command from the Terminal.  Currently, `NGStoolkit` supports RNA-seq analysis via a single-click automated RNA-seq pipeline starting with `bz2` files through differential expression.

There are 6 input parameters:

*	Pipeline type { `simple`, `complete` }.  `simple` is for individual flow cells where all people want are the stats for a given run, and `complete` is through from `bz2` to differential expression.
*	Genome type { `hg19`, `hg38`, `mm10`, `rn6` }.  Specify the genome for the samples.
*	Read type { `paired`, `single` }.  Paired/single ended read flag.
*	Configuration file { `<filename>` }.  This file is a list of sample IDs and the groups for differential expression.  First column is the sample ID and the second column is the group.
*	Stranded information { `u`, `s`, `r` }.  Type of strandedness to take the counts from STAR output.  `u` = unstranded, `s` = stranded, `r` = reverse stranded as defined by STAR.

## Usage

`sh rna_pipeline.sh <pipeline_type> <genome> <paired/single> <config> <strandedness>`

## Future Plans

Current work is underway to expand `NGStoolkit` into ChIP-seq, methyl-seq, and exome-seq analysis territory.  We are also expanding from a hardcoded HPC environment to a flexible cloud infrastructure.

## Funding

`NGStoolkit` is an ongoing bioinformatics software project financially supported by the
United States Department of Defense (DoD) through the National Defense Science and Engineering
Graduate Fellowship (NDSEG) Program. This research was conducted with Government support under
and awarded by DoD, Army Research Office (ARO), National Defense Science and Engineering
Graduate (NDSEG) Fellowship, 32 CFR 168a.
