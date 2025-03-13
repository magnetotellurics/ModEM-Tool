%
%  script for compiling readedi & writeedi.
%
mex -fortran module_PLTIO_v1.1.F90 readedi.F90 -output readedi
mex -fortran module_PLTIO_v1.1.F90 writeedi.F90 -output writeedi
