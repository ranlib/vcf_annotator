#!/usr/bin/env python3
"""
Purpose:
get annotation for vcf file from EnsEmbl via rest calls

Usage:
usage: vcf_annotator_post.py [-h] -i VCF_INPUT -o TSV_OUTPUT -j JSON_OUTPUT [-r SERVER] [-v] [-q] [-m]

get annotation for vcf file

optional arguments:
  -h, --help            show this help message and exit
  -i VCF_INPUT, --input VCF_INPUT
                        input vcf file
  -o TSV_OUTPUT, --output TSV_OUTPUT
                        output tsv file
  -j JSON_OUTPUT, --json JSON_OUTPUT
                        output json file
  -r SERVER, --server SERVER
                        EnsEMBL server to use (default=https://grch37.rest.ensembl.org)
  -v, --verbose         turn on debugging output
  -q, --query           turn on rest calls to server for gene information
  -m, --maf             turn on rest calls to server for MAF information

Comments:
1) -q parameter: turn on at first run in order to produce the json file with the annotation information
   for a vcf file with about 11,000 rows this takes about 20 min with bulk rest call, 200 entries at a time
2) -m parameter: turn on if Minor Allele Frequency should be added to output tsv file,
   rest call for one entry at a time, for vcf file with about 11,000 rows this takes about 2 hours
   one-at-a-time calls because no all variants have annotation information

Notes:
1. Depth of sequence coverage at the site of variation.
2. Number of reads supporting the variant.
3. Percentage of reads supporting the variant versus those supporting reference reads.
4. Using the VEP hgvs API, get the
   a) gene of the variant,
   b) type of variation (substitution, insertion, CNV, etc.) and
   c) their effect (missense, silent, intergenic, etc.).
   The API documentation is available at: https://rest.ensembl.org/#VEP
5. The minor allele frequency of the variant if available.
6. Any additional annotations that you feel might be relevant.

Tools used:
https://vcfpy.readthedocs.io/en/stable/index.html
https://pypi.org/project/nested-lookup/

server = "https://grch37.rest.ensembl.org"

Author: Dieter Best
Date: 03. September 2022

"""
import os
import sys
import argparse
import logging
import json
import itertools
import time
import csv

import vcfpy
import requests
from nested_lookup import nested_lookup


def fetch_endpoint(server: str, request: str):
    """
    fetch endpoint data
    return json
    """
    try:
        response = requests.get(
            server + request, headers={"Accept": "application/json"}
        )
    except requests.HTTPError as exception:
        print(exception)
        return None

    if not response.ok:
        logging.warning("response from server not ok")
        return None

    return response.json()


def fetch_endpoint_post(
    server: str,
    hgvs_list,
    ext: str = "/vep/human/hgvs"
):
    """
    fetch all data from endpoint
    return json
    """
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    data = {}
    data["hgvs_notations"] = hgvs_list
    try:
        response = requests.post(server + ext, headers=headers, json=data)
    except requests.HTTPError as exception:
        print(exception)
        return None

    if not response.ok:
        logging.warning("response from server not ok")

    return response.json() if response.ok else None


def get_chunk(alist: [], size_of_chunk: int):
    """
    divide list into list of sublists
    each sublist has length size_of_chunk
    """
    for i in range(0, len(alist), size_of_chunk):
        yield alist[i : i + size_of_chunk]


def vcf_annotator_post(params):
    """
    get annotation for variant via REST from EnsEmbl

    Input:
    vcf_input: str: input vcf file
    verbose: bool: turn on verbose debugging output

    Output:
    tsv_output: str: output tsv file
    """

    # setup vcf reader
    reader = vcfpy.Reader.from_path(params.vcf_input)

    # loop over variants
    # fill list of hgvs records for rest call
    hgvs_list = []  # list of hgvs records
    variants = {}  # list of variants, each entry dictionary for variant information
    for record in reader:
        # get information from record needed to make rest call
        alt_allele_list = [allele.value for allele in record.ALT]
        number_of_reads_with_variant_list = record.INFO["TR"]

        # loop over alt_alleles
        for allele, number_of_reads in zip(
            alt_allele_list, number_of_reads_with_variant_list
        ):
            hgvs = f"{record.CHROM}:g.{record.POS}{record.REF}>{allele}"
            hgvs_list.append(hgvs)

            # storage for variant information
            variants_dict = {
                "chromosome": record.CHROM,
                "position": record.POS,
                "ref_allele": record.REF,
                "alt_allele": allele,
                "total_number_of_reads": record.INFO["TC"],
                "number_of_reads_with_variant": number_of_reads,
                "percent_reads_with_variant": round(
                    float(number_of_reads) * 100 / float(record.INFO["TC"]), 2
                ),
            }
            variants[hgvs] = variants_dict

    # REST calls to EnsEmbl to get the data for variant
    # limitation of 200 entries at a time
    if params.query:
        json_list = []
        counter = 0
        size = 200  # upper limit for EnsEmbl, if bigger then error
        for chunk in list(get_chunk(hgvs_list, size)):
            if params.verbose:
                print(counter)
            start = time.time()
            data = fetch_endpoint_post(params.server, chunk)
            end = time.time()
            print("Time elapsed:", end - start)
            if params.verbose:
                print(json.dumps(data, indent=4))
            if data is not None:
                json_list.extend(data)
            counter += 1

        # write data to json file
        with open(params.json_output, "w") as json_file:
            json.dump(json_list, json_file, indent=4)
    else:  # open json file
        with open(params.json_output, "r") as json_file:
            json_list = json.load(json_file)

    # process variant data
    # add information to variants list
    # write annotation tsv file
    keys = [
        "gene_symbol",
        "impact",
        "most_severe_consequence",
        "id",
        "consequence_terms",
        "MAF",
    ]
    for hgvs in variants:
        for key in keys:
            variants[hgvs][key] = "Null"

    # Note: not all entries in variants have annotation information
    for entry in json_list:
        hgvs = nested_lookup("input", entry)[0]
        for key in ["gene_symbol", "impact", "most_severe_consequence", "id"]:
            variants[hgvs][key] = ",".join(set(nested_lookup(key, entry)))

        variants[hgvs]["consequence_terms"] = ",".join(
            set(itertools.chain(*nested_lookup("consequence_terms", entry)))
        )

    # MAF
    # get rsids
    # get json data with MAF from EnsEmbl
    # some entries have rsid = Null
    if params.maf:
        json_population_list = []
        counter = 0
        for hgvs in variants:
            if variants[hgvs]["id"] != "Null":
                counter += 1
                rsid = [id for id in variants[hgvs]["id"].split(",") if "rs" in id]
                if len(rsid) > 0:
                    if params.verbose:
                        print(counter)
                    start = time.time()
                    request = f"/variation/human/{rsid[0]}?pop=1"
                    population = fetch_endpoint(params.server, request)
                    end = time.time()
                    print("Time elapsed:", end - start)
                    if params.verbose:
                        print(json.dumps(population, indent=4))
                    if population is not None:
                        json_population_list.append(population)
                    maf = nested_lookup("MAF", population)
                    if len(maf) > 0:
                        if maf[0] is not None:
                            variants[hgvs]["MAF"] = maf[0]
            else:
                logging.warning("no rsid for %s", hgvs)
        # write data to json file
        with open(params.json_population_output, "w") as json_file:
            json.dump(json_population_list, json_file, indent=4)

    # write variants to tsv file
    with open(params.tsv_output, "w") as csv_file:
        fields = variants[next(iter(variants))].keys()
        writer = csv.DictWriter(csv_file, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        for key, value in variants.items():
            writer.writerow(value)

    return 0


if __name__ == "__main__":
    # get input vcf file
    DESCRIPTION = "get annotation for vcf file"
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument(
        "-i",
        "--input",
        dest="vcf_input",
        type=str,
        help="input vcf file",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="tsv_output",
        type=str,
        help="output tsv file",
        required=True,
    )
    parser.add_argument(
        "-j",
        "--json",
        dest="json_output",
        type=str,
        help="output json file",
        required=True,
    )
    parser.add_argument(
        "-p",
        "--population",
        dest="json_population_output",
        type=str,
        help="output json file for population information",
        required=True,
    )
    parser.add_argument(
        "-r",
        "--server",
        dest="server",
        type=str,
        help="EnsEMBL server to use (default=%(default)s)",
        default="https://grch37.rest.ensembl.org",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        help="turn on debugging output",
    )
    parser.add_argument(
        "-q",
        "--query",
        dest="query",
        action="store_true",
        help="turn on rest calls to server for gene information",
    )
    parser.add_argument(
        "-m",
        "--maf",
        dest="maf",
        action="store_true",
        help="turn on rest calls to server for MAF information",
    )
    args = parser.parse_args()

    # check input vcf file exists
    if not os.path.exists(args.vcf_input):
        logging.error("input file %s does not exist.", args.vcf_input)
        sys.exit(1)

    # check if json file exists
    if args.query is None:
        if not os.path.exists(args.json_output):
            logging.error("input file %s does not exist, run with -q option.", args.json_output)
            sys.exit(1)

    ERR = vcf_annotator_post(args)

    sys.exit(ERR)
