%
%  Script for testing edmod3d.m.
%  Written by Bo Yang, IGG, CUG, and CEOAS, OSU, 2013/11/04.
%  Last change: 2013/11/07 10:53:00.
%
clc
clear
m2d.xn = -4:.2:4;
m2d.zn = -2:.2:2;
[xx,zz] = meshgrid(m2d.xn,m2d.zn);
v = xx.*exp(-xx.^2-zz.^2);
m2d.pp  = v;
out = edmod2d(m2d);
