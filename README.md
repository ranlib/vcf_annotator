# vcf_annotator
annotate vcf file via rest calls to EMBL

Purpose:
get annotation for vcf file from EnsEmbl via rest calls

```
Usage:
vcf_annotator_post.py [-h] -i VCF_INPUT -o TSV_OUTPUT -j JSON_OUTPUT [-r SERVER] [-v] [-q] [-m]

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
```

Documentation:
1) -q parameter: turn on at first run in order to produce the json
   file with the annotation information for a vcf file with about
   11,000 rows this takes about 20 min with bulk rest call, 200
   entries at a time
2) -m parameter: turn on if Minor Allele Frequency should be added to
   output tsv file, rest call for one entry at a time, for vcf file
   with about 11,000 rows this takes about 2 hours one-at-a-time calls
   because no all variants have annotation information
3) data structure that hold the information and information is added
   to variants = dictionary index on hgvs terms, each entry =
   dictionary of annotation information
4) procedure has 6 steps:

   1) read in vcf file
   2) add information from vcf file to variants data structure
   3) get via REST calls to EnsEmbl annotation information in json
   format for hgvs terms in variants data structure
   4) parse json and add information to variants data strucutre
   5) get rsid information from variants data structure and via REST call to EnsEmbl get MAF
   6) write variants data structure (dictionary of dictionaries) to tsv file

Notes:
1. Depth of sequence coverage at the site of variation.
2. Number of reads supporting the variant.
3. Percentage of reads supporting the variant versus those supporting reference reads.
4. Using the VEP hgvs API, get the
   1) gene of the variant,
   2) type of variation (substitution, insertion, CNV, etc.) and
   3) their effect (missense, silent, intergenic, etc.).
   The API documentation is available at: https://rest.ensembl.org/#VEP
5. The minor allele frequency of the variant if available.
6. Any additional annotations that you feel might be relevant.

Tools used:
1) https://vcfpy.readthedocs.io/en/stable/index.html
2) https://pypi.org/project/nested-lookup/

server = "https://grch37.rest.ensembl.org"

