# Usage
# -----
# docker build -t rnapromo:dev -f run-rnapromo.dockerfile .


FROM ubuntu:16.04
ENV DEBIAN_FRONTEND=noninteractive


RUN apt -y update && apt -y install \
    build-essential \
    clang \
    imagemagick \
    libstdc++5 \
    make \
    unzip \
    zip \
    && apt clean all


# Install old perl
WORKDIR /
ADD https://github.com/Perl/perl5/archive/refs/tags/v5.20.3-RC2.tar.gz /
RUN tar -xzf v5.20.3-RC2.tar.gz \
    && cd perl5-5.20.3-RC2 \
    && ./Configure -des -Dprefix=/usr \
    && make -j \
    && make install


# Install xerces
ADD https://dlcdn.apache.org//xerces/c/3/sources/xerces-c-3.2.5.tar.gz /
RUN tar -xzf xerces-c-3.2.5.tar.gz \
    && cd xerces-c-3.2.5 \
    && ./configure \
    && make -j \
    && make install


# Install RNApromo
WORKDIR /
ADD rnapromo.zip .
RUN unzip rnapromo.zip


# Entrypoint
WORKDIR /rnapromo
ENTRYPOINT ["perl"]
CMD ["--version"]
