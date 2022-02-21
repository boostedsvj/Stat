#!/bin/bash

./impacts.sh 3100 /store/user/pedrok/SVJ2017/Limits/datacards_Sep13/ "-t -1 --expectSignal 0" asimovSig0
./impacts.sh 3100 /store/user/pedrok/SVJ2017/Limits/datacards_Sep13/ "-t -1 --expectSignal 1.6" asimovSigRexp cut
./impacts.sh 3100 /store/user/pedrok/SVJ2017/Limits/datacards_Sep13/ "-t -1 --expectSignal 0.2" asimovSigRexp bdt

