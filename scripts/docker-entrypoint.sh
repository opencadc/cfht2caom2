#!/bin/bash

if [[ ! -e ${PWD}/cache.yml ]]
then
  cp /cache.yml ${PWD}
fi

if [[ ! -e ${PWD}/config.yml ]]
then
  cp /config.yml ${PWD}
fi

if [[ ! -e ${PWD}/state.yml ]]; then
  echo "${@}"
  if [[ "${@}" == "cfht_run_state" ]]; then
    echo "if statement true"
    yesterday=$(date -d yesterday "+%d-%b-%Y %H:%M")
    echo "bookmarks:
    cfht_timestamp:
      last_record: $yesterday
" > ${PWD}/state.yml
  fi
fi

exec "${@}"
