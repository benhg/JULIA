# JULIA

## What is this?

JULIA is a project designed to explain/verify taxon source within a large and complex dataset of possibly contaminated transcriptomes. By "contaminated," here we mean we are not sure exactly from where each of our transcripts came.

## How does it work?

JULIA uses Parsl as an execution environment to carry out alignment and indexing at scale.

## How do I use it?

1. Locate a database of raw reads
2. Create a FASTA file of input transcripts
3. Call JULIA using python and Parsl
4. Extract results with JULIA's data-extraction application.

## What does it do?

Using Alex Dobin's STAR to perform alignment, JULIA indexes and aligns transcripts against a database of raw reads. From there, it is able to generate information about taxon source for each of the transcripts input to it.

## Who wrote it?

Ellen Richards and Ben Glick, expanding on an earlier project by Julia Somers, with advisory support from Greta Binford.

## Who is it named after?

While JULIA is an acronym, it is additionally named after Julia Somers, Lewis & Clark Class of 2019, who designed the project that preceeded it.
