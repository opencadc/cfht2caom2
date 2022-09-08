FROM debian:buster-slim
ENV DEBIAN_FRONTEND=noninteractive

ADD https://www.python.org/ftp/python/3.9.10/Python-3.9.10.tgz /usr/local/src/

RUN apt-get update --no-install-recommends \
    && apt-get install -y \
        gcc \
        g++ \
        git \
        libc6-dev \
        libcfitsio-bin \
        libgdbm-dev \
        libncursesw5-dev \
        libreadline-gplv2-dev \
        libsqlite3-dev \
        libssl-dev \
        libbz2-dev \
        libffi-dev \
        libtool \
        make \
        saods9 \
        tk-dev \
        xvfb \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/ /tmp/* /var/tmp/*

RUN cd /usr/local/src \
    && tar zxvf Python-3.9.10.tgz \
    && cd Python-3.9.10 \
    && ./configure --enable-optimizations --prefix=/usr/local \
    && make \
    && make install \
    && ln -s /usr/local/bin/python3 /usr/local/bin/python \
    && ln -s /usr/local/bin/pip3 /usr/local/bin/pip

RUN pip install --no-cache-dir wheel

RUN pip install --no-cache-dir "astropy<5" \
    && pip install pytz \
    && pip install pyyaml

WORKDIR /usr/src/app

ARG OPENCADC_BRANCH=master
ARG OPENCADC_REPO=opencadc

RUN git clone https://github.com/opencadc/caom2tools.git && \
    cd caom2tools && \
    pip install ./caom2utils && \
    cd ..

RUN pip install git+https://github.com/${OPENCADC_REPO}/caom2pipe@${OPENCADC_BRANCH}#egg=caom2pipe

RUN pip install git+https://github.com/${OPENCADC_REPO}/cfht2caom2@${OPENCADC_BRANCH}#egg=cfht2caom2

RUN useradd --create-home --shell /bin/bash cadcops
RUN chown -R cadcops:cadcops /usr/src/app
USER cadcops

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]
