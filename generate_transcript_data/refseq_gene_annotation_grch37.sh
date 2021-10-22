#!/bin/bash

pyreference_args=()

# Having troubles with corrupted files downloading via FTP from NCBI via IPv6, http works ok

filename=ref_GRCh37.p5_top_level.gff3.gz
url=http://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/BUILD.37.3/GFF/${filename}
pyreference_file=${filename}.json.gz
if [[ ! -e ${filename} ]]; then
  wget ${url}
fi
if [[ ! -e ${pyreference_file} ]]; then
  pyreference_gff_to_json.py --url "${url}" --gff3 "${filename}"
fi
pyreference_args+=(--pyreference-json ${pyreference_file})


filename=ref_GRCh37.p9_top_level.gff3.gz
url=http://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.103/GFF/${filename}
pyreference_file=${filename}.json.gz
if [[ ! -e ${filename} ]]; then
  wget ${url}
fi
if [[ ! -e ${pyreference_file} ]]; then
  pyreference_gff_to_json.py --url "${url}" --gff3 "${filename}"
fi
pyreference_args+=(--pyreference-json ${pyreference_file})


filename=ref_GRCh37.p10_top_level.gff3.gz
url=http://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.104/GFF/${filename}
pyreference_file=${filename}.json.gz
if [[ ! -e ${filename} ]]; then
  wget ${url}
fi
if [[ ! -e ${pyreference_file} ]]; then
  pyreference_gff_to_json.py --url "${url}" --gff3 "${filename}"
fi
pyreference_args+=(--pyreference-json ${pyreference_file})


filename=ref_GRCh37.p13_top_level.gff3.gz
url=http://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.105/GFF/${filename}
pyreference_file=${filename}.json.gz
if [[ ! -e ${filename} ]]; then
  wget ${url}
fi
if [[ ! -e ${pyreference_file} ]]; then
  pyreference_gff_to_json.py --url "${url}" --gff3 "${filename}"
fi
pyreference_args+=(--pyreference-json ${pyreference_file})


# These all have the same name, so rename them based on release ID
for release in 105.20190906 105.20201022; do
  filename=GCF_000001405.25_GRCh37.p13_genomic.${release}.gff.gz
  url=http://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/annotation_releases/${release}/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gff.gz
  pyreference_file=${filename}.json.gz
  if [[ ! -e ${filename} ]]; then
    wget ${url} --output-document=${filename}
  fi
  if [[ ! -e ${pyreference_file} ]]; then
    pyreference_gff_to_json.py --url "${url}" --gff3 "${filename}"
  fi
  pyreference_args+=(--pyreference-json ${pyreference_file})
done


merged_file="pyhgvs_transcripts_refseq_grch37.json.gz"
if [[ ! -e ${merged_file} ]]; then
  BASE_DIR=$(dirname ${BASH_SOURCE[0]})

  python3 ${BASE_DIR}/pyreference_to_pyhgvs_json.py ${pyreference_args[@]} --output ${merged_file}
fi
