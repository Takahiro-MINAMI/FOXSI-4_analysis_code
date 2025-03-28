#!/bin/bash

rm en.txt
touch en.txt
root -l -b -q peaktotxt_am.cpp
root -l -b -q peaktotxt_ba.cpp
root -l -b -q peaktotxt_fe.cpp
root -l -b -q pedestal_peaktotxt.cpp
root -l -b -q genCalTree.cpp

open enCalDataTree.pdf
