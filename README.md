# FastMatrix #

FastMatrix is an exercise from a college Linear Algebra class, circa 2007.

The FastMatrix class is a math-oriented template C++ class for storing arbitrary sized matrices of arbitrary numeric type (64-bit double precision floating point
numbers are used by default to represent real entries). fastMAtrix may also be instantiated with other types such as complex numbers. The goal of the fastMatrix class is to offer a number of common mathematical operations in a
reasonably efficient manner.

Copyright 2016 Matthew Bennett, original code 2007

Licensed under MPL 2.0. See LICENSE for details.

## Usage Instructions ##

Library does multiplication, transpose, compute inverse, DFT, FFT, QR factorization, LU decomposition, extract eigenvectors, etc. 

Please Read fastMatrix.pdf for details.

## Compiling ##
All you have to do is compile the cpp files using g++ or any ansi standard compiler.

example: 

g++ -o fastMatrix *.cpp *.h

Type ./fastMatrix to run it. You may wish to modify main.cpp to experiment with the various features of the fastMatrix class. The interface is documented in fastMatrix.h and within the manual.pdf.

