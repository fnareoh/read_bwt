Read-optimized Burrow Wheeler transform
=========

###Brief description
-----------
This is a implementation of a Read-optmized BWT a joint work of Travis Gagie, Garance Gourdel, Gonzalo Navarro, Nicola Prezza, Giovanna Rosone, Jared Simpson.
The implementation was done by Garance Gourdel (garance.gourdel@ens-paris-saclay.fr).

###Requirements
------------
The use of the data stucture requires [the sdsl-lite library][https://github.com/simongog/sdsl-lite].

###Acknowlegment
------------
All the file in the `/internal` directory are the data structure Nicola Prezza implemented for the [R-index][https://github.com/nicolaprezza/r-index]. The `data/dna.in` file used by `extract.cpp` to generate the reads is the dna file from [Pizza & Chili][http://pizzachili.dcc.uchile.cl/].

###Installation
------------
```sh
git clone git@github.com:fnareoh/read_bwt.git
cd read_bwt
```
Then you should modify the sdsl path in the `Makefile`.
```sh
make all
```

###Getting Started
------------
If you want to test the data structure defined in `wt.cpp`, `wt.hpp`, `read_bwt.cpp` and `read_bwt.hpp` you can just `make all` and run :
```sh
./extract && ./extension && ./new_extension && ./test_struct
```
If you wish to use the data structures you should copy the files `wt.cpp`, `wt.hpp`, `read_bwt.cpp` and `read_bwt.hpp` and then use the structures

read_bwt is initialized with the concatenation of all the read and potential extension and a bit-vector marking the position of the extensions or dolars.

Warning: There is a problem with the structure if you want to use the structure, you should initialize the structure save it and load it again before using it.
