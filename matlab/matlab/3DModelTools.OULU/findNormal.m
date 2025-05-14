function [A,X] = findNormal(A,B)
% Finds normal vector AX of the line AB (at the point A)
%
% INPUT
% A, B 2d vectors
%
% OUTPUT
% A, X 2d vectors

% 1st: Move the origin to the point A
B(1) = B(1)-A(1);
B(2) = B(2)-A(2);

% 2nd: Rotate +90 degrees
Rot = [cos(pi/2) -sin(pi/2); sin(pi/2) cos(pi/2)];
N = Rot*B';

% 3rd: Move back to the original coordinates
X(1) = N(1)+A(1);
X(2) = N(2)+A(2);