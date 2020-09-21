#! /usr/bin/env bash

cd ../build

# kill all jobs if ctr-c is received
trap 'kill %1; kill %2; kill %3; kill %4; kill %5; kill %6' SIGINT

./run_norkyst &
./run_norkyst_ref &
./run_nordic &
./run_nordic_ref &
./run_arctic &
./run_arctic_ref &

wait
