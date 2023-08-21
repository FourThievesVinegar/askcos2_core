FROM continuumio/miniconda3:23.3.1-0

RUN apt update && apt -y install gcc g++ make

COPY env.yaml /tmp/env.yaml

RUN conda env update -n base -f /tmp/env.yaml

COPY . /askcos2_core
WORKDIR /askcos2_core
