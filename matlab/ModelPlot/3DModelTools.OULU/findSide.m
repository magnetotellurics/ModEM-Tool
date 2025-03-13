function res = findSide(Q,A,B)
% Finds out which side of the line AB is the point Q.
%
% INPUT
% Q, A, B 2d vectors
%
% OUTPUT
% -1 = "over" the line
%  0 = on the line
%  1 = "under" the line

% Does this work?
area = (B(1)-A(1))*(Q(2)-A(2))-(B(2)-A(2))*(Q(1)-A(1));
res = sign(area);