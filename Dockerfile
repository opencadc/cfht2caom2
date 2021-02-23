FROM opencadc/matplotlib:3.8-slim
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y && apt-get dist-upgrade -y && \
    apt-get install -y \
    xvfb \
    git \
    python3-astropy \
    python3-pip \
    python3-tz \
    python3-yaml \
    saods9 && \
    rm -rf /var/lib/apt/lists/ /tmp/* /var/tmp/*

RUN pip install  --no-cache-dir \
        aplpy \
        bs4 \
        cadcdata \
        cadctap \
        caom2 \
        caom2repo \
        caom2utils \
        ftputil \
        importlib-metadata \
        PyYAML \
        pytz \
        spherical-geometry \
        vos

RUN rm -rf /var/lib/apt/lists/ /tmp/* /var/tmp/*

WORKDIR /usr/src/app

RUN git clone https://github.com/HEASARC/cfitsio && \
  cd cfitsio && \
  ./configure --prefix=/usr && \
  make -j 2 && \
  make shared && \
  make install && \
  make fitscopy && \
  cp fitscopy /usr/local/bin && \
  make clean && \
  cd /usr/src/app

ARG OPENCADC_BRANCH=master
ARG OPENCADC_REPO=opencadc

RUN git clone --depth=1 https://github.com/${OPENCADC_REPO}/caom2pipe.git --branch=${OPENCADC_BRANCH} && \
    rm -rf ./caom2pipe/.git && \
    pip install ./caom2pipe

RUN git clone --depth=1 https://github.com/${OPENCADC_REPO}/cfht2caom2.git --branch=${OPENCADC_BRANCH} && \
    rm -rf ./cfht2caom2/.git && \
    cp ./cfht2caom2/scripts/config.yml / && \
    cp ./cfht2caom2/scripts/cache.yml / && \
    cp ./cfht2caom2/scripts/docker-entrypoint.sh / && \
    pip install ./cfht2caom2

ENTRYPOINT ["/docker-entrypoint.sh"]

