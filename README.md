# fDE

## About
fortran differential evolution and other parallel samplers

Algorithm https://link.springer.com/article/10.1007/s11222-008-9104-9

Implemented with reference to python - spotpy - demcz and R - BayesianTools - mcmcDEzs 

Motivation : reimplemented in Fortran for use in CARDAMOM 

## Obtaining

### Binaries
TODO

### Build from source
git clone

```
cd fDE
mkdir build
cd build
cmake ..
make 
```

#### Developer

```
cmake .. -DENABLE_TESTING
make
ctest .
```

## Using
It's a library.  It has no main exectuable.  

It exposes functions:

### Example
```fortran
use fDE
```
