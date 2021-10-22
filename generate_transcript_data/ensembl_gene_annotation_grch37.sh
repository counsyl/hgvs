#!/bin/bash

# v81 (points to 75) and earlier at GTFs that don't have transcript versions - just skip them

#82 is first GFF3 for GRCh37
#83 has no data
#84 is 82 again
#86 is 85 again
pyreference_args=()
for release in 82 85 87; do
  filename=Homo_sapiens.GRCh37.${release}.gff3.gz
  url=ftp://ftp.ensembl.org/pub/grch37/release-${release}/gff3/homo_sapiens/${filename}
  pyreference_file=${filename}.json.gz
  if [[ ! -e ${filename} ]]; then
    wget ${url}
  fi
  if [[ ! -e ${pyreference_file} ]]; then
    pyreference_gff_to_json.py --url "${url}" --gff3 "${filename}"
  fi
  pyreference_args+=(--pyreference-json ${pyreference_file})
done

merged_file="pyhgvs_transcripts_ensembl_grch37.json.gz"
if [[ ! -e ${merged_file} ]]; then
  BASE_DIR=$(dirname ${BASH_SOURCE[0]})

  python3 ${BASE_DIR}/pyreference_to_pyhgvs_json.py ${pyreference_args[@]} --output ${merged_file}
fi
