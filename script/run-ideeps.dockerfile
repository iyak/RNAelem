# Usage
# -----
# git clone https://github.com/xypan1232/iDeepS.git
# docker build -t ideeps:dev -f run-ideeps.dockerfile iDeepS


FROM ubuntu:16.04

WORKDIR /iDeepS
COPY * .
RUN chmod a+x RNAshapes
ENV PATH="${PATH}:/iDeepS"

RUN apt-get update && apt-get install -y ca-certificates
RUN update-ca-certificates
RUN apt-get --allow-insecure-repositories update; exit 0
RUN apt-get --allow-unauthenticated install --yes\
    build-essential\
    cmake\
    curl\
    gfortran\
    git-core\
    git\
    graphviz-dev\
    libblas-dev\
    libbz2-dev\
    libc-dev\
    libcairo2-dev\
    libeigen3-dev\
    libffi-dev\
    libfreetype6-dev\
    libjpeg62-dev\
    liblapack-dev\
    liblzma-dev\
    libopenbabel-dev\
    libpng-dev\
    openbabel\
    pkg-config \
    python2.7-dev\
    software-properties-common\
    tar\
    unzip\
    vim\
    wget\
    zip\
    &&\
    apt-get remove -y --purge libzmq-dev software-properties-common libc-dev libreadline-dev && \
    apt-get autoremove -y && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN wget https://sourceforge.net/projects/swig/files/swig/swig-2.0.0/swig-2.0.0.tar.gz
RUN tar -xvf swig-2.0.0.tar.gz
WORKDIR /iDeepS/swig-2.0.0
RUN ./configure && make -j && make install

WORKDIR /iDeepS
RUN curl https://bootstrap.pypa.io/pip/2.7/get-pip.py --output get-pip.py
RUN python get-pip.py
RUN pip install --upgrade setuptools pip
RUN pip install numpy==1.6.0\
    && pip install matplotlib==1.3.1\
    && pip install --upgrade numpy==1.15
RUN pip install scipy==0.14.0
RUN pip install\
    biopython==1.66\
    cvxopt==1.1.7\
    cvxpy==0.2.12\
    dill==0.2.1\
    ecos==1.0.3\
    esmre==0.3.1\
    flake8==2.4.1\
    git+https://github.com/nils-werner/pymf.git@9d3b8d49941b018473de01cb28b97d42b8675857\
    joblib==0.8.4\
    keras==1.1.2\
    mpld3==0.2\
    networkx==1.9\
    nose==1.3.7\
    pandas==0.14.1\
    pybedtools==0.7.4\
    pygraphviz==1.2\
    regex==2014.01.10\
    reportlab==3.2.0\
    requests==2.7.0\
    scikit-learn==0.17.0\
    scikit-neuralnetwork==0.3\
    scripttest==1.3\
    scs==1.0.1\
    seaborn==0.5\
    tensorflow==1.14.0\
    theano==0.9.0\
    weblogo==3.4

RUN unzip EDeN.zip
WORKDIR /iDeepS/EDeN
RUN pip install -e . --no-deps

WORKDIR /iDeepS
ENTRYPOINT ["python", "ideeps.py"]
