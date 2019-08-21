FROM nfcore/base
LABEL author="onur.yukselen@umassmed.edu" description="Docker image containing all requirements for the dolphinnext/chipseq pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
# Install standard utilities for HOMER
RUN apt-get -y install zip unzip gcc g++ make
# Install dolphin-tools
RUN git clone https://github.com/UMMS-Biocore/tools /usr/local/bin/dolphin-tools
RUN mkdir -p /project /nl /mnt /share
ENV PATH /opt/conda/envs/dolphinnext-chipseq-1.0/bin:/usr/local/bin/dolphin-tools/:$PATH
