#!/bin/sh

# Requires Intel FORTRAN and C compilers for parrallel and O3 optimizations.
# (segfaults produced by gfortran and gcc).

icc -c pipelibc_g77.c
ifort -c pipelibfort_g77.f
ifort -O3 -parallel ptrinum.f pipelibfort_g77.o pipelibc_g77.o -o ptrinum
ifort -O3 -parallel ptrigather.f pipelibfort_g77.o pipelibc_g77.o -o ptrigather
ifort -O3 -parallel ptrilink.f pipelibfort_g77.o pipelibc_g77.o -o ptrilink
