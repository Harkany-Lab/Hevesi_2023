# Start from nvidia-docker image with drivers pre-installed to use a GPU
FROM nvcr.io/nvidia/pytorch:22.12-py3

LABEL maintainer="Evgenii O. Tretiakov <evgenii.tretiakov@meduniwien.ac.at>"

# Install curl and sudo and git and miniconda and pytorch, cudatoolkit, pytables, and cellbender
RUN apt-get update && apt-get install -y --no-install-recommends \
	curl \
	ca-certificates \
	sudo \
	&& sudo apt-get install -y --no-install-recommends \
	git \
	bzip2 \
	&& sudo rm -rf /var/lib/apt/lists/*

ENV PATH=/opt/conda/bin:$PATH
ENV CONDA_AUTO_UPDATE_CONDA=false

ENV DOCKER='true'

RUN git clone https://github.com/broadinstitute/CellBender.git cellbender \
	&& yes | python -m pip install -e cellbender \
	&& sudo rm -rf ~/.cache/pip

