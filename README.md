# cardamom-samplers

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
cd cardamom-samplers
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
It's a library.  It has no main exectuable other than the test/ examples.

Use the samplers in a fortran program.

It exposes functions ``MHMCMC``, ``DEMCz``:

### Example
