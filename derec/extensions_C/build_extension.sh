#!/bin/bash

gcc lxnorm.c -Wall -fpic -O2 -c -o lxnorm.o
gcc -shared -o liblxnorm.so lxnorm.o 
