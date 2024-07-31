# Usage
# -----
# docker build -t graphprot2:dev -f run-graphprot2.dockerfile .


FROM pytorch/pytorch:1.2-cuda10.0-cudnn7-devel
ENV FORCE_CUDA=1


RUN curl https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1604/x86_64/3bf863cc.pub \
    | apt-key add -
RUN apt update && apt install -y \
    build-essential \
    git \
    vim \
    && rm -rf /var/lib/apt/lists/*


RUN pip install --no-cache --no-cache-dir --force-reinstall \
    --find-links https://pytorch-geometric.com/whl/torch-1.5.0.html \
    logomaker==0.8 \
    markdown==3.2.2 \
    torch-cluster==1.4.4 \
    torch-geometric==1.3.1 \
    torch-scatter==1.3.1 \
    torch-sparse==0.4.0 \
    torch-spline-conv==1.1.0 \
    ushuffle


RUN conda install --channel bioconda \
    bedtools \
    seaborn==0.10.1 \
    ucsc-bigwigaverageoverbed \
    ucsc-twobitinfo \
    ucsc-twobittofa \
    viennarna=2.4.14


# Install GraphProt2
RUN git clone --branch master --single-branch https://github.com/BackofenLab/GraphProt2.git \
    && cd GraphProt2 \
    && git checkout 9e02890abafe6c506c64b8552b4540927e5ce0c0 \
    && python -m pip install . --ignore-installed --no-deps -vv


# Entrypoint
ENTRYPOINT ["graphprot2"]
CMD ["--help"]

