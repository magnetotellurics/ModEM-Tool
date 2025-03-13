function [mOut] = cTimesMod(c,mIn)
%  Initial implementation of scalar multiplication of a model
%   space object by a real scalar; no error type checking;
%   scalar has to be first argument
%
%  Usage: [mOut] = cTimesMod(c,mIn)

mOut = mIn;
mOut.v = c*mOut.v;
