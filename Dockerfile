FROM opencadc/matplotlib:3.9-slim
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update --no-install-recommends && apt-get dist-upgrade -y && \
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
        importlib-metadata \
        python-dateutil \
        PyYAML \
        spherical-geometry \
        vos

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
ARG PIPE_BRANCH=master
ARG PIPE_REPO=opencadc

RUN git clone https://github.com/opencadc/cadctools.git && \
    cd cadctools && \
    pip install ./cadcutils && \
    pip install ./cadcdata && \
    cd ..

RUN pip install git+https://github.com/${OPENCADC_REPO}/caom2pipe@${OPENCADC_BRANCH}#egg=caom2pipe

RUN pip install git+https://github.com/${PIPE_REPO}/cfht2caom2@${PIPE_BRANCH}#egg=cfht2caom2

WORKDIR /usr/src/app

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]

