FROM ubuntu:latest
RUN apt-get update --yes
RUN apt-get install --yes python3 python3-pip
RUN pip3 install vcfpy pysam nested_lookup requests
COPY ./vcf_annotator_post.py /usr/local/bin



