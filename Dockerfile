FROM opencadc/astropy:4-3.7-alpine
  
RUN apk --no-cache add \
        bash \
        coreutils \
        gcc \
        g++ \
        libmagic

RUN pip install cadcdata && \
        pip install cadctap && \
        pip install caom2 && \
        pip install caom2repo && \
        pip install caom2utils && \
        pip install deprecated && \
        pip install ftputil && \
        pip install importlib-metadata && \
        pip install pytz && \
        pip install PyYAML && \
        pip install spherical-geometry && \
        pip install vos

WORKDIR /usr/src/app

RUN pip install bs4

RUN git clone https://github.com/SharonGoliath/caom2tools.git && \
    pip install ./caom2tools/caom2 && \
    pip install ./caom2tools/caom2utils

RUN git clone https://github.com/SharonGoliath/caom2pipe.git && \
    pip install ./caom2pipe

COPY ./cfht2caom2_unit/scripts/docker-entrypoint.sh /usr/src/app/cfht2caom2
COPY ./cfht2caom2_unit/scripts/config.yml /config.yml
COPY ./cfht2caom2_unit/scripts/cache.yml /cache.yml

ENTRYPOINT ["/usr/src/app/cfht2caom2/docker-entrypoint.sh"]

