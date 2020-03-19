FROM debian:buster-slim
# ADD docker-apt-install /usr/local/bin
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y && apt-get dist-upgrade -y

RUN apt-get install -y \
    xvfb \
    git \
    python3-astropy \
    python3-pip \
    python3-tz \
    python3-yaml

RUN pip3 install  --no-cache-dir \
        cadcdata \
        cadctap \
        caom2 \
        caom2repo \
        caom2utils \
        deprecated \
        ftputil \
        importlib-metadata \
        spherical-geometry \
        vos

RUN apt-get install -y saods9

RUN rm -rf /var/lib/apt/lists/ /tmp/* /var/tmp/*

RUN git clone https://github.com/HEASARC/cfitsio && \
  cd cfitsio && \
  ./configure --prefix=/usr && \
  make -j 2 && \
  make shared && \
  make install && \
  make fitscopy && \
  cp fitscopy /usr/local/bin && \
  make clean

WORKDIR /usr/src/app

RUN pip3 install bs4

ARG OPENCADC_REPO=opencadc
ARG OMC_REPO=opencadc-metadata-curation

RUN git clone https://github.com/${OPENCADC_REPO}/caom2tools.git && \
    pip3 install ./caom2tools/caom2 && \
    pip3 install ./caom2tools/caom2utils

RUN git clone https://github.com/${OMC_REPO}/caom2pipe.git && \
    pip3 install ./caom2pipe

RUN git clone https://github.com/${OMC_REPO}/cfht2caom2.git && \
    cp ./cfht2caom2/scripts/config.yml / && \
    cp ./cfht2caom2/scripts/docker-entrypoint.sh / && \
    pip3 install ./cfht2caom2

ENTRYPOINT ["/docker-entrypoint.sh"]

