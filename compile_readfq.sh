#!/bin/bash
gcc -std=c99 -Wall -Doff64_t=__off64_t -lz readfq.c -o readfq
