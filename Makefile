CC=g++
CFLAGS= -std=c++11
SDSLFLAGS1= -DNDEBUG -O3  -I/home/fnareoh/include -L/home/fnareoh/lib
SDSLFLAGS0= -g -O0  -I/home/fnareoh/include -L/home/fnareoh/lib
SDSLFLAGS2=  -lsdsl -ldivsufsort -ldivsufsort64
.DEFAULT_GOAL := all

extract: extract.cpp
	$(CC) -o extract extract.cpp

program: program.cpp
	$(CC) $(CFLAGS) $(SDSLFLAGS0) program.cpp  -o program $(SDSLFLAGS2)

test: test.cpp
	$(CC) $(CFLAGS) $(SDSLFLAGS0) test.cpp -o test_program $(SDSLFLAGS2)

test_struct: test_struct.cpp
	$(CC) $(CFLAGS) $(SDSLFLAGS0) test_struct.cpp -o test_struct $(SDSLFLAGS2)

extension: extension.cpp
	$(CC) $(CFLAGS) $(SDSLFLAGS1) extension.cpp -o extension $(SDSLFLAGS2)

new_extension: new_extension.cpp
		$(CC) $(CFLAGS) $(SDSLFLAGS1) new_extension.cpp -o new_extension $(SDSLFLAGS2)

all:
	make extract extension new_extension test_struct

clean:
	rm extract extension new_extension test_struct
	rm data/*.out data/reads* data/pattern.in
