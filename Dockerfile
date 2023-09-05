FROM mambaorg/micromamba:1.4.7

USER root
# Keep the base environment activated
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN apt update && apt -y install gcc g++ make

# Use micromamba to resolve conda-forge, much faster than conda
RUN micromamba install -y python=3.10.12 pip=23.2.1 rdkit=2023.03.3 -c conda-forge
RUN pip install \
    celery==5.2.7 \
    fastapi==0.95.1 \
    gevent==22.10.2 \
    pandas==1.5.3 \
    "passlib[bcrypt]"==1.7.4 \
    pillow==9.5.0 \
    protobuf==3.19.0 \
    pydantic==1.10.12 \
    pymongo==4.4.1 \
    pytest==7.4.1 \
    "python-jose[cryptography]"==3.3.0 \
    python-multipart==0.0.6 \
    pyyaml==6.0.1 \
    redis==4.3.6 \
    requests==2.31.0 \
    svgutils==0.3.4 \
    tqdm==4.66.1 \
    uvicorn==0.21.1

COPY . /ASKCOSv2/askcos2_core
WORKDIR /ASKCOSv2/askcos2_core
