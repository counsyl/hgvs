"""
    Converts PyReference JSON.gz files (created from RefSeq/Ensembl GTF or GFFs) into format easy to load w/PyHGVS

    @see https://bitbucket.org/sacgf/pyreference/

    Dave Lawrence (davmlaw@gmail.com) on 22/10/2021
"""
import gzip
import json
from argparse import ArgumentParser
from typing import Dict


def handle_args():
    parser = ArgumentParser(description='Convert multiple PyReference json.gz files into one for PyHGVS')
    parser.add_argument('--pyreference-json', required=True, nargs="+", action="extend",
                        help='PyReference JSON.gz - list OLDEST to NEWEST (newest is kept)')
    parser.add_argument('--output', required=True, help='Output filename')
    return parser.parse_args()


def convert_gene_pyreference_to_gene_version_data(gene_data: Dict) -> Dict:
    gene_version_data = {
        'biotype': ",".join(gene_data["biotype"]),
        'description': gene_data.get("description"),
        'gene_symbol': gene_data["name"],
    }

    if hgnc_str := gene_data.get("HGNC"):
        # Has HGNC: (5 characters) at start of it
        gene_version_data["hgnc"] = hgnc_str[5:]

    # Only Ensembl Genes have versions
    if version := gene_data.get("version"):
        gene_data["version"] = version

    return gene_version_data


def convert_transcript_pyreference_to_pyhgvs(transcript_data: Dict) -> Dict:
    start = transcript_data["start"]
    end = transcript_data["stop"]
    strand = transcript_data["strand"]
    # PyHGVS has cds_start/cds_end be equal to start/end for non-coding transcripts
    cds_start = transcript_data.get("cds_start", start)
    cds_end = transcript_data.get("cds_end", end)
    # PyHGVS exons are in genomic order, PyReference are in stranded
    features = transcript_data["features_by_type"]
    exons = [[ed["start"], ed["stop"]] for ed in features["exon"]]
    cdna_match = []
    for cdm in features.get("cDNA_match", []):
        if "cdna_strand" in cdm:  # None in Human RefSeq so haven't handled
            raise ValueError("Haven't handled stranded Target alignment")
        cdna_match.append([cdm.get(k) for k in ["start", "stop", "cdna_start", "cdna_end", "gap"]])

    if strand == '-':
        exons.reverse()
        cdna_match.reverse()

    pyhgvs_data = {
        'chrom': transcript_data["chr"],
        'start': start,
        'end': end,
        'strand': strand,
        'cds_start': cds_start,
        'cds_end': cds_end,
        'exons': exons,
    }

    # Optional stuff
    if cdna_match:
        pyhgvs_data["cdna_match"] = cdna_match
    if transcript_data.get("partial"):
        pyhgvs_data["partial"] = 1
    other_chroms = transcript_data.get("other_chroms")
    if other_chroms:
        pyhgvs_data["other_chroms"] = other_chroms

    # Few extra fields (pyHGVS doesn't currently use)
    pyhgvs_data["biotype"] = ",".join(transcript_data["biotype"])
    return pyhgvs_data


def main():
    args = handle_args()

    gene_versions = {}  # We only keep those that are in the latest transcript version
    transcript_versions = {}

    for pyref_filename in args.pyreference_json:
        print(f"Loading '{pyref_filename}'")
        with gzip.open(pyref_filename) as f:
            pyref_data = json.load(f)

            url = pyref_data["reference_gtf"]["url"]

            # PyReference stores transcripts under genes, while PyReference only has transcripts (that contain genes)
            transcript_gene_version = {}

            for gene_id, gene in pyref_data["genes_by_id"].items():
                if version := gene.get("version"):
                    gene_accession = f"{gene_id}.{version}"
                else:
                    gene_accession = gene_id

                gene_version = convert_gene_pyreference_to_gene_version_data(gene)
                gene_version["url"] = url
                gene_versions[gene_accession] = gene_version

                for transcript_accession in gene["transcripts"]:
                    transcript_gene_version[transcript_accession] = gene_accession

            for transcript_accession, transcript in pyref_data["transcripts_by_id"].items():
                hgvs_data = convert_transcript_pyreference_to_pyhgvs(transcript)
                gene_accession = transcript_gene_version[transcript_accession]
                gene_version = gene_versions[gene_accession]
                hgvs_data["id"] = transcript_accession
                hgvs_data["gene_version"] = gene_accession
                hgvs_data["gene_name"] = gene_version["gene_symbol"]
                hgvs_data["url"] = url
                transcript_versions[transcript_accession] = hgvs_data

    print("Writing PyHGVS data")
    with gzip.open(args.output, 'w') as outfile:
        data = {
            "transcripts": transcript_versions,
            "genes": gene_versions,
        }
        json_str = json.dumps(data)
        outfile.write(json_str.encode('ascii'))


if __name__ == '__main__':
    main()
