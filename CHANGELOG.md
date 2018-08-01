# HGVS library change log

## 0.9.5 (2018-07-01)
  -  Fix pip installs by not using the internal API
## 0.9.4  / 0.9.3 (2015-03-11)
  - `NC_000005.10:g.177421339_177421327delACTCGAGTGCTCC` appears in ClinVar, and is an invalid name (the genomic start/stop coords are not in increasing order). This causes parse_hgvs_name to raise an IndexError. It should raise InvalidHGVSName instead
## 0.9.2 (2014-12-11)
  - Rename package to pyhgvs to avoid naming conflict with other libraries.

## 0.9.1 (2014-03-10)
  - Normalize variants by default when generating names.

## 0.9 (2014-01-16)
  - Improved testing.
