#!/bin/bash
gcc -std=c99 -Wall -Doff64_t=__off64_t -lz kliass.c -o kliass
