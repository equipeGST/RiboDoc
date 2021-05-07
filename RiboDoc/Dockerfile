FROM ubuntu:18.04

LABEL version="0.6"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y wget git-all \
    && rm -rf /var/lib/apt/lists/*

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm Miniconda3-latest-Linux-x86_64.sh

RUN git clone https://github.com/equipeGST/RiboDoc.git

ENV BASH_ENV ~/.bashrc

ENV PATH=/root/miniconda3/bin:$PATH
ENV PATH=/root/miniconda3/envs/RiboDoc_env/bin:$PATH

SHELL ["/bin/bash", "-c"]

RUN conda env create -f /RiboDoc/RiboDoc/RiboDoc_env.yml \
    && conda clean -a -y
RUN conda init bash
RUN echo "conda activate RiboDoc_env" >> /root/.bashrc \
    && conda info --envs \
    && conda list --name RiboDoc_env

RUN R -e "library('devtools');devtools::install_github('LabTranslationalArchitectomics/riboWaltz', dependencies = FALSE);"

CMD ["bash","/RiboDoc/RiboDoc/RiboDoc.sh"]
