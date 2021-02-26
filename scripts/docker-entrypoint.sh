#!/bin/bash

if [[ ! -e ${PWD}/cache.yml ]]
then
  cp /usr/local/bin/cache.yml ${PWD}
fi

if [[ ! -e ${PWD}/config.yml ]]
then
  cp /usr/local/bin/config.yml ${PWD}
fi

if [[ ! -e ${PWD}/state.yml ]]; then
  if [[ "${@}" == "cfht_run_state" ]]; then
    yesterday=$(date -d yesterday "+%d-%b-%Y %H:%M")
    echo "bookmarks:
    cfht_timestamp:
      last_record: $yesterday
" > ${PWD}/state.yml
  fi
fi

exec "${@}"
