FROM continuumio/miniconda3

RUN conda install numpy scipy>=1 networkx>=2 pymongo psutil coverage \
    && conda install -c conda-forge spglib \
    && pip install periodictable seekpath \ 
    && conda clean --all

RUN apt-get update \
    && apt-get install -y openmpi-bin
