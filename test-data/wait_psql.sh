#!/bin/bash

for ((i=0;i<100;i++))
do
    DB_CONNECTABLE=$(nc -z localhost 5432 >/dev/null 2>&1; echo "$?")
        if [[ $DB_CONNECTABLE -eq 0 ]]; then
                break
        fi
    sleep 3
done

if ! [[ $DB_CONNECTABLE -eq 0 ]]; then
        echo "Cannot connect to database"
    exit "${DB_CONNECTABLE}"
fi
