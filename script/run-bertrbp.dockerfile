# Usage
# -----
# docker build -t bert-rbp:dev -f run-bertrbp.dockerfile .


FROM pytorch/pytorch:1.5-cuda10.1-cudnn7-runtime


RUN apt update && apt install -y \
    build-essential \
    git \
    libbz2-dev \
    make \
    vim \
    zip \
    && rm -rf /var/lib/apt/lists/*


RUN pip install --no-cache \
    absl-py==2.1.0 \
    asn1crypto==1.3.0 \
    backcall==0.1.0 \
    beautifulsoup4==4.9.0 \
    biopython==1.81 \
    boto3==1.33.13 \
    botocore==1.33.13 \
    cachetools==5.3.3 \
    certifi==2020.4.5.1 \
    cffi==1.14.0 \
    chardet==3.0.4 \
    click==8.1.7 \
    conda-build==3.18.11 \
    conda-package-handling==1.6.0 \
    conda==4.8.3 \
    cryptography==2.8 \
    cycler==0.11.0 \
    decorator==4.4.2 \
    filelock==3.0.12 \
    fonttools==4.38.0 \
    fsspec==2023.1.0 \
    gdown==4.7.3 \
    glob2==0.7 \
    google-auth-oauthlib==0.4.6 \
    google-auth==2.32.0 \
    grpcio==1.62.2 \
    huggingface-hub==0.16.4 \
    idna==2.8 \
    importlib-metadata==6.7.0 \
    ipython-genutils==0.2.0 \
    ipython==7.13.0 \
    jedi==0.16.0 \
    Jinja2==2.11.1 \
    jmespath==1.0.1 \
    joblib==1.3.2 \
    kiwisolver==1.4.5 \
    libarchive-c==2.9 \
    lief==0.9.0 \
    Markdown==3.4.4 \
    MarkupSafe==2.1.5 \
    matplotlib==3.5.3 \
    mkl-fft==1.0.15 \
    mkl-random==1.1.0 \
    mkl-service==2.3.0 \
    numpy==1.18.1 \
    oauthlib==3.2.2 \
    olefile==0.46 \
    packaging==24.0 \
    pandas==1.3.5 \
    parso==0.6.2 \
    patsy==0.5.6 \
    pexpect==4.8.0 \
    pickleshare==0.7.5 \
    Pillow==7.0.0 \
    pkginfo==1.5.0.1 \
    prompt-toolkit==3.0.4 \
    protobuf==3.20.0 \
    psutil==5.7.0 \
    ptyprocess==0.6.0 \
    pyahocorasick==2.0.0 \
    pyasn1-modules==0.3.0 \
    pyasn1==0.5.1 \
    pybedtools==0.10.0 \
    pycosat==0.6.3 \
    pycparser==2.19 \
    Pygments==2.6.1 \
    pyOpenSSL==19.1.0 \
    pyparsing==3.1.2 \
    pysam==0.22.1 \
    PySocks==1.7.1 \
    python-dateutil==2.9.0.post0 \
    pytz==2019.3 \
    PyYAML==5.3.1 \
    regex==2024.4.16 \
    requests-oauthlib==2.0.0 \
    requests==2.22.0 \
    rsa==4.9 \
    ruamel-yaml==0.15.87 \
    s3transfer==0.8.2 \
    sacremoses==0.0.53 \
    scikit-learn==0.22.2 \
    scipy==1.7.3 \
    seaborn==0.12.2 \
    sentencepiece==0.1.91 \
    seqeval==1.2.2 \
    six==1.14.0 \
    soupsieve==2.0 \
    statsmodels==0.13.5 \
    tensorboard-data-server==0.6.1 \
    tensorboard-plugin-wit==1.8.1 \
    tensorboard==2.11.2 \
    tensorboardX==2.6.2.2 \
    tokenizers==0.19.1 \
    torch==1.5.0 \
    torchvision==0.6.0a0+82fd1c8 \
    tqdm==4.42.1 \
    traitlets==4.3.3 \
    typing-extensions==4.7.1 \
    urllib3==1.25.8 \
    wcwidth==0.1.9 \
    Werkzeug==2.2.3 \
    zipp==3.15.0


# Install BERT-RBP
WORKDIR /root/github/kkyamada/
RUN git clone --branch master --single-branch \
    https://github.com/kkyamada/bert-rbp.git \
    && cd bert-rbp \
    && git checkout abf624d5f60247a9e4921fc110bff93e441fc37a \
    && gdown -O 3-new-12w-0.zip 1nVBaIoiJpnwQxiz4dSq6Sv9kBKfXhZuM \ 
    && unzip -o 3-new-12w-0.zip


# Entrypoint
WORKDIR /root/github/kkyamada/bert-rbp/scripts
ENTRYPOINT ["bash", "train_and_test.sh"]
CMD ["DATA", "../3-new-12w-0"]
