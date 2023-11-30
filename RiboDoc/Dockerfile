FROM ubuntu:20.04

LABEL version="0.9.0"

ENV LC_ALL=C

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y wget git-all \
    && rm -rf /var/lib/apt/lists/*

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py310_23.5.2-0-Linux-x86_64.sh \
    && bash Miniconda3-py310_23.5.2-0-Linux-x86_64.sh -b -u -p /bin/miniconda3/ \
    && rm Miniconda3-py310_23.5.2-0-Linux-x86_64.sh

RUN git clone https://github.com/equipeGST/RiboDoc.git

RUN chmod +rx /RiboDoc/RiboDoc/RiboDoc.sh

RUN mv /root/.bashrc /bin/ && ln -s /bin/.bashrc /root/.bashrc && \
    mv /root/.conda /bin/ && ln -s /bin/.conda /root/.conda && \
    mv /root/.profile /bin/ && ln -s /bin/.profile /root/.profile && \
    ln -s /root/ /bin/

ENV PATH=/bin/miniconda3/bin:$PATH
ENV PATH=/bin/miniconda3/envs:$PATH
ENV PATH=$PATH:/bin

SHELL ["/bin/bash", "-c"]

RUN conda env create -f /RiboDoc/RiboDoc/RiboDoc_env.yml
RUN conda init bash

RUN echo "conda activate RiboDoc_env" >> /bin/.bashrc

ENV PATH=/bin/miniconda3/envs/RiboDoc_env/bin:$PATH

RUN /bin/miniconda3/envs/RiboDoc_env/bin/R -e "library('devtools');devtools::install_github('LabTranslationalArchitectomics/riboWaltz@v1.2.0', dependencies = FALSE);"

WORKDIR /data/

ENTRYPOINT ["bash","/RiboDoc/RiboDoc/RiboDoc.sh"]
