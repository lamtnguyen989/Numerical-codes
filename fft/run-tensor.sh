#!/bin/bash

mkdir -p bin/

make tensor
./bin/tensor-fft > run-tensor.log
cat run-tensor.log