#!/bin/bash

set -e

# VERSION=`cat VERSION.txt`

docker build -t brownmp/virtect:devel .
