#!/bin/bash
for i in {1..50}
do
    echo ./mce97110h.exe
    ./mce97110h.exe
    echo h2root mchc12ntp.hbook $i.root
    h2root mchc12ntp.hbook $i.root
done
