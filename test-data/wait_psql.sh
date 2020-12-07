#!/bin/bash

# Consider there's only one docker running
DOCKER_ID=`docker ps --format={{.Names}} | head -n1`

echo "Watching docker id $DOCKER_ID"

for ((i=0;i<100;i++))
do
    DB_CONNECTABLE=$(docker logs $DOCKER_ID | grep -q 02-search.sh >/dev/null 2>&1; echo "$?")
        if [[ $DB_CONNECTABLE -eq 0 ]]; then
                break
        fi
    sleep 3
done

sleep 10

if ! [[ $DB_CONNECTABLE -eq 0 ]]; then
        echo "Cannot connect to database"
        docker logs $DOCKER_ID
    exit "${DB_CONNECTABLE}"
fi
