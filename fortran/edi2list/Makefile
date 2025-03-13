# Makefile for compiling and testing edi2list.
# Written by Bo Yang, 2014/03/20.
#
FC=gfortran
FLAG=-O2
BIN=edi2list
SRC=edi2list_main.f90
MOD=module_LATLON_v1.0.f90 \
	module_PLTIO_v1.1.f90 \
	ParseCode.f90

all: $(BIN) tst

$(BIN): $(SRC) $(MOD)
	$(FC) $(FLAG) -o $(BIN) $(MOD) $(SRC)

tst:
	cp $(BIN) ./test
	cd ./test;./$(BIN) runfile

clean:
	-rm *.mod *.o
