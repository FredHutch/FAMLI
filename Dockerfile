FROM ubuntu:16.04
MAINTAINER sminot@fredhutch.org

# Install prerequisites
RUN apt update && \
	apt-get install -y build-essential wget unzip python2.7 \
					   python-dev git python-pip bats awscli curl \
					   libcurl4-openssl-dev make gcc zlib1g-dev

# Set the default langage to C
ENV LC_ALL C

# Use /share as the working directory
RUN mkdir /share
WORKDIR /share

# Add files
RUN mkdir /usr/famli
ADD requirements.txt /usr/famli

# Install python requirements
RUN pip install -r /usr/famli/requirements.txt && rm /usr/famli/requirements.txt


# Install DIAMOND v0.9.10
RUN cd /usr/famli && \
	wget -q https://github.com/bbuchfink/diamond/releases/download/v0.9.10/diamond-linux64.tar.gz && \
	tar xzf diamond-linux64.tar.gz && \
	mv diamond /usr/bin/ && \
	rm diamond-linux64.tar.gz


# Install the SRA toolkit
RUN cd /usr/local/bin && \
	wget -q https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.2/sratoolkit.2.8.2-ubuntu64.tar.gz && \
	tar xzf sratoolkit.2.8.2-ubuntu64.tar.gz && \
	ln -s /usr/local/bin/sratoolkit.2.8.2-ubuntu64/bin/* /usr/local/bin/ && \
	rm sratoolkit.2.8.2-ubuntu64.tar.gz


# Make a user to install Aspera Connect
RUN useradd -ms /bin/bash user && \
	usermod -g root user && \
	chmod g+w /usr/local/bin
USER user

# Install Aspera Connect
RUN cd /usr/local/bin && \
	wget -q https://download.asperasoft.com/download/sw/connect/3.7.4/aspera-connect-3.7.4.147727-linux-64.tar.gz && \
	tar xzvf aspera-connect-3.7.4.147727-linux-64.tar.gz && \
	bash aspera-connect-3.7.4.147727-linux-64.sh

# Set the ASPERA_KEY environment variable
ENV ASPERA_KEY=/home/user/.aspera/connect/etc/asperaweb_id_dsa.openssh

USER root
RUN ln -s /home/user/.aspera/connect/bin/* /usr/local/bin/

# Install the FASTX Toolkit
RUN cd /usr/local/bin && \
	wget -q http://hannonlab.cshl.edu/fastx_toolkit/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2 && \
	tar xf fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2 && \
	rm fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64.tar.bz2 && \
	mv bin/* ./


# Add the run script to the PATH
ADD famli.py /usr/famli
ADD famli /usr/famli/famli
RUN cd /usr/famli && \
	ln -s /usr/famli/famli.py /usr/bin/
ENV PYTHONPATH="/usr/famli:${PYTHONPATH}"


# Run tests and then remove the folder
ADD tests /usr/famli/tests
RUN bats /usr/famli/tests/ && rm -r /usr/famli/tests/
