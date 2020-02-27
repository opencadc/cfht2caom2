# FROM opencadc/astropy:3.6-alpine

FROM python:3.7

# RUN apk --no-cache add \
#     bash \
#     coreutils \
#     git
RUN pip install numpy && \
    pip install astropy

RUN pip install cadcdata && \
    pip install cadctap && \
    pip install caom2 && \
    pip install caom2repo && \
    pip install caom2utils && \
    pip install deprecated && \
    pip install ftputils && \
    pip install PyYAML && \
    pip install pytz && \
    pip install spherical-geometry && \
    pip install vos

# RUN apk --no-cache add libxt-dev
# RUN pip install pyds9

RUN wget http://ds9.si.edu/download/debian10/ds9.debian10.8.1.tar.gz && \
    tar -xvf ds9.debian10.8.1.tar.gz && \
    cp ds9 /usr/local/bin

WORKDIR /usr/src/app

# RUN git clone https://github.com/opencadc-metadata-curation/caom2pipe.git && \
#     pip install ./caom2pipe
  
RUN git clone https://github.com/opencadc-metadata-curation/cfht2caom2.git && \
  cp ./cfht2caom2/scripts/config.yml / && \
  cp ./cfht2caom2/scripts/cache.yml / && \
  cp ./cfht2caom2/scripts/docker-entrypoint.sh / && \
  pip install ./cfht2caom2

ENTRYPOINT ["/docker-entrypoint.sh"]

