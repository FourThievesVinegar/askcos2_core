FROM continuumio/miniconda3:22.11.1-alpine
COPY env.yaml /tmp/env.yaml

RUN conda env create -f /tmp/env.yaml
RUN echo "conda activate askcos2-core" >> ~/.profile

COPY . /askcos2_core
WORKDIR /askcos2_core
