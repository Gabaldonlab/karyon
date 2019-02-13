#!/bin/bash
###
# Karyon Test .
# version 0.1a
###

current=$(pwd)
docker stop karyon_test
docker run -dit --name=karyon_test -v $current:/home/ --rm ubuntu
docker exec -it -w /home karyon_test sh /home/install.sh
docker exec -it -w /home karyon_test /bin/bash

echo "Ended!"