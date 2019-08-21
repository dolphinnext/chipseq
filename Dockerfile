FROM nfcore/base
LABEL author="onur.yukselen@umassmed.edu" description="Docker image containing all requirements for the dolphinnext/chipseq pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
# Install standard utilities for HOMERv4.10
RUN apt-get clean all
RUN apt-get update
RUN apt-get dist-upgrade -y
RUN apt-get -y install zip unzip gcc g++ make
# Install dolphin-tools
RUN git clone https://github.com/dolphinnext/tools /usr/local/bin/dolphin-tools
RUN mkdir -p /project /nl /mnt /share
ENV PATH /opt/conda/envs/dolphinnext-chipseq-1.0/bin:/usr/local/bin/dolphin-tools/:$PATH
