function [Ax] = multA(x)

NORMALIZE = 1;
y = Jfreq(x);
Ax = JTfreq(y,NORMALIZE);
Ax = Ax(1,:);
