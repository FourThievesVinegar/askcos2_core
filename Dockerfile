FROM continuumio/miniconda3:22.11.1-alpine

COPY env.yaml /tmp/env.yaml

RUN conda env update -n base -f /tmp/env.yaml

COPY . /askcos2_core
WORKDIR /askcos2_core
