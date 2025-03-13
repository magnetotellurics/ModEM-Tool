Feb 20, 2014
1, run "mex -setup" to setup the fortran compiler for your matlab;
2, run script "compile_writereadedi.m" to compile "readedi.f90 & writeedi.f90", you will get "readedi.mex & writeedi.mex";
3, run edirazor -> click "Load list" -> open "testedi.list" to load an EDI file list with the full path of each edi file;
4, choose the site in the list box on the right of the axis;
5, use left btn to remove period at the right of the cursor, use the right btn to bring it back;
6, click "Save as..." btn to save the selection info in a new edi file for each site. There will be one new keyword called ">FSEL"(short for frequency selection) append to the EDI files.


Tips:
1, If the mex compiler complain about some strange fortran format errors, it might need to remove
the "/fixed" option in "mexopts.bat" file. "mexopts.bat" located at "C:\Users\Administrator\AppData\Roaming\MathWorks\MATLAB\R2011a\mexopts.bat" on Windows 7 x64 system.