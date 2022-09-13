FROM ubuntu:18.04

LABEL version="0.8.3"

ENV LC_ALL=C

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y wget git-all \
    && rm -rf /var/lib/apt/lists/*

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -u -p /bin/miniconda3/ \
    && rm Miniconda3-latest-Linux-x86_64.sh

RUN git clone https://github.com/equipeGST/RiboDoc.git

RUN chmod +rx /RiboDoc/RiboDoc/RiboDoc.sh

RUN mv /root/.bashrc /bin/ && ln -s /bin/.bashrc /root/.bashrc && \
    mv /root/.conda /bin/ && ln -s /bin/.conda /root/.conda && \
    mv /root/.profile /bin/ && ln -s /bin/.profile /root/.profile && \
    ln -s /root/ /bin/

ENV PATH=/bin/miniconda3/bin:$PATH
ENV PATH=/bin/miniconda3/envs/RiboDoc_env/bin:$PATH
ENV PATH=$PATH:/bin

SHELL ["/bin/bash", "-c"]

RUN conda env create -f /RiboDoc/RiboDoc/RiboDoc_env.yml
RUN conda init bash
RUN echo "conda activate RiboDoc_env" >> /bin/.bashrc \
    && conda info --envs \
    && conda list --name RiboDoc_env

RUN R -e "library('devtools');devtools::install_github('LabTranslationalArchitectomics/riboWaltz@v1.2.0', dependencies = FALSE);"

ENTRYPOINT ["bash","/RiboDoc/RiboDoc/RiboDoc.sh"]
