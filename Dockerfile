FROM opencadc/astropy:3.6-alpine

RUN apk --no-cache add \
    bash \
    coreutils \
    git

RUN pip install cadcdata && \
    pip install cadctap && \
    pip install caom2 && \
    pip install caom2repo && \
    pip install caom2utils && \
    pip install PyYAML && \
    pip install spherical-geometry && \
    pip install vos

WORKDIR /usr/src/app

RUN git clone https://github.com/opencadc-metadata-curation/caom2pipe.git && \
    pip install ./caom2pipe
  
RUN git clone https://github.com/opencadc-metadata-curation/cfht2caom2.git && \
  cp ./cfht2caom2/scripts/config.yml / && \
  cp ./cfht2caom2/scripts/docker-entrypoint.sh / && \
  pip install ./cfht2caom2

RUN apk --no-cache del git

ENTRYPOINT ["/docker-entrypoint.sh"]

