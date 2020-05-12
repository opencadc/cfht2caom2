#!/bin/bash

# always over-write the operational cache.yml file, since there's no
# other way to make sure version-controlled additions are made operational
cp /cache.yml ${PWD}

if [[ ! -e ${PWD}/config.yml ]]
then
  cp /config.yml ${PWD}
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
