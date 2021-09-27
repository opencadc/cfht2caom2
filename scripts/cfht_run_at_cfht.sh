#!/bin/bash

COLLECTION="cfht"
IMAGE="cfht_build"

# echo "Get a proxy certificate"
# cp $HOME/.ssl/data/makanaproxy.pem ./ || exit $?

echo "Run image ${IMAGE}"
docker run --init --rm --mount "type=bind,src=${PWD},dst=/usr/src/app" --mount "type=bind,src=/data/makana,dst=/data" ${IMAGE} ${COLLECTION}_run_state

date
exit 0
