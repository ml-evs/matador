FROM continuumio/miniconda3

RUN conda install numpy scipy pymongo psutil coverage \
    && conda install -c conda-forge spglib=1.9.9.18 \
    && pip install periodictable \ 
    && conda clean --all
