#!/bin/bash

#run on small vcf file for testing
#time vcf_annotator_post.py -i small_test_vcf_data.vcf -o output_post_small_test_vcf_data.tsv -j output_post_small_test_vcf_data.json -p output_post_small_test_vcf_population.json -v -q -m 

# run with docker
time docker run --rm -v $PWD:/mnt -w /mnt dbest/vcf_annotator:v1.0 vcf_annotator_post.py -i small_test_vcf_data.vcf -o output_post_small_test_vcf_data.tsv -j output_post_small_test_vcf_data.json -p output_post_small_test_vcf_population.json -v -q -m

# the big file
#time vcf_annotator_post.py -i test_vcf_data.txt -o output_post_test_vcf_data.tsv -j output_post_test_vcf_data.json -p output_post_test_vcf_population.json -v -m
