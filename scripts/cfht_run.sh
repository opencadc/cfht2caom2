#!/bin/bash

COLLECTION="cfht"
IMAGE="opencadc/${COLLECTION}2caom2"

echo "Get a proxy certificate"
cp $HOME/.ssl/cadcproxy.pem ./ || exit $?

echo "Get image ${IMAGE}"
docker pull ${IMAGE}

echo "Run image ${IMAGE}"
docker run --init --rm -v ${PWD}:/usr/src/app/ --user $(id -u):$(id -g) -e HOME=/usr/src/app ${IMAGE} ${COLLECTION}_run || exit $?

date
exit 0
