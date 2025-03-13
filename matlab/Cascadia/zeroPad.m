function [cPad] = zeroPad(c)

[n,m] = size(c);
c = [ c zeros(n,1)];
cPad = [ c ; zeros(1,m+1)];
