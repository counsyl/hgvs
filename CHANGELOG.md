# HGVS library change log

## Unreleased
  - Use same code for calculating transcript position of start/stop codons.
  - Be able to pass in pre-calculated start/stop from transcript_json if available

## 0.12.1 (2021-12-09)
  - Fix issue #61 - Regex for reference HGVS without reference base

## 0.12.0 (2021-11-24)
  - cDNA gaps
  - Refactor models
  - Support non-coding, mitochondria and LRG transcripts
  - Validate coordinate span equals reference length
  - Scripts to generate transcripts from RefSeq and Ensembl GTFs on the web

## 0.9.2 (2014-12-11)
  - Rename package to pyhgvs to avoid naming conflict with other libraries.

## 0.9.1 (2014-03-10)
  - Normalize variants by default when generating names.

## 0.9 (2014-01-16)
  - Improved testing.
