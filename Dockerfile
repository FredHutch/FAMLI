FROM ubuntu:24.04
LABEL MAINTAINER=sminot@fredhutch.org

# Install prerequisites
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
	apt-get install -y build-essential wget unzip python3.12 \
					   python3-dev git python3-pip bats \
					   python3.12-venv \
					   libcurl4-openssl-dev make gcc zlib1g-dev curl

# Set the default language to C
ENV LC_ALL=C

# Use /share as the working directory
RUN mkdir /share
WORKDIR /share

# Add /scratch
RUN mkdir /scratch

# Folder for installation
RUN mkdir /usr/famli

# Install DIAMOND v0.9.10
RUN cd /usr/local/bin && \
	wget -q https://github.com/bbuchfink/diamond/releases/download/v0.9.22/diamond-linux64.tar.gz && \
	tar xzf diamond-linux64.tar.gz && \
	rm diamond-linux64.tar.gz


# Install the SRA toolkit
RUN cd /usr/local/bin && \
	wget -q https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.2/sratoolkit.2.8.2-ubuntu64.tar.gz && \
	tar xzf sratoolkit.2.8.2-ubuntu64.tar.gz && \
	ln -s /usr/local/bin/sratoolkit.2.8.2-ubuntu64/bin/* /usr/local/bin/ && \
	rm sratoolkit.2.8.2-ubuntu64.tar.gz


# Install the FASTX Toolkit
RUN cd /usr/local/bin && \
	wget -q http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2 && \
	tar xf fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2 && \
	rm fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2 && \
	mv bin/* ./

# Make a virtual environment for Python 3.12
RUN python3.12 -m venv /usr/envs/famli
RUN /bin/bash -c "source /usr/envs/famli/bin/activate"
ENV PATH="/usr/envs/famli/bin:$PATH"

# Install UV
RUN pip install uv

# Install FAMLI from PyPI
RUN uv pip install famli==1.2.2

# Add the local directory to the container
ADD . /usr/famli

# Add the taxonomic analysis script to the PATH
RUN ln -s /usr/famli/diamond-tax.py /usr/local/bin/

# Run tests and then remove the folder
RUN bats /usr/famli/tests/ && rm -r /usr/famli/tests/
