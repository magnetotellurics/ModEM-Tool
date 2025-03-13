1, run "mex -setup" to setup the fortran compiler for your matlab;
2, run script "compile.m" to compile "dellatlon.f90", you will get "dellatlon.mex";
3, use "dellatlon" as a regular matlab function. type help dellatlon for usage information.

Tips:
1, If the mex compiler complain about some strange fortran format errors, it might need to remove
the "/fixed" option in "mexopts.bat" file. "mexopts.bat" located at "C:\Users\Administrator\AppData\Roaming\MathWorks\MATLAB\R2011a\mexopts.bat" on Windows 7 x64 system.