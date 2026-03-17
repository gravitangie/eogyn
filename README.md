# eogyn: <ins>*e*</ins>ffective-<ins>*o*</ins>ne-body <ins>*g*</ins>eneric d<ins>*yn*</ins>amics 

[![DOI](https://zenodo.org/badge/1184559551.svg)](https://doi.org/10.5281/zenodo.19073497) 

Welcome! Here you can find `eogyn` - a C code that solves the dynamics of black hole binaries with precessing spins on generic orbits within the effective-one-body (EOB) formalism. Our work is based on the Hamiltonian presented in [this paper](https://inspirehep.net/literature/1395081).

## Getting started

**Prerequisites**: C compiler

Clone the repository:

```
git clone https://github.com/gravitangie/eogyn.git
cd eogyn
```

Compile the code:

```
cd C
make
```

## Usage

Run the code with an input parfile:

```
./eogyn -p parfiles/example.par
```

## Developers

This code was developed by Angelica Albertini and Michal Stratený under the supervision of Georgios Loukes-Gerakopoulos.

## Acknowledgements

This code was developed within the project No. 107324 funded by the Grant Agency of Charles University.