#!/bin/bash

mkdir -p bin/

make basic
./bin/fft > run.log
cat run.log