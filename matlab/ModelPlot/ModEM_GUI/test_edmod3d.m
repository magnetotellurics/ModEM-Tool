%
%  Script for testing edmod3d.m.
%  Written by Bo Yang, IGG, CUG, and CEOAS, OSU, 2013/11/04.
%  Last change: 2013/11/05 20:24:25.
%%
clc
clear
m3d.xn = -2:.2:2;
m3d.yn = -3:.2:3;
m3d.zn = -4:.2:4;
[xx,yy,zz] = meshgrid(m3d.xn,m3d.yn,m3d.zn);
v = xx.*exp(-xx.^2-yy.^2-zz.^2);
m3d.v  = v;
outmod = edmod3d(m3d);