#!/bin/bash

export name_list="cli engine_ser mdworker"
export pid_list;

function get_pidlist()
{
    pid_list=`ps | grep -E 'cli|engine_ser|mdworker' | grep -v grep | awk '{print $1}'`
}

get_pidlist
while [[ ${pid_list} ]]; do
    killall ${name_list}
    get_pidlist
done
echo "Killed visit processes ${name_list}"
