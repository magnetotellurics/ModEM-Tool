function res = GCD(A,B,varargin)
% The Great-circle distance of points A(Lat,Long) and B(Lat,Long)
% with radius 1 or (varargin). 
% note that lenght(A) = lenght(B) or length(B) = 1.
%
% Numerically (more) stabile Haversine formula.

if nargin > 2
    r = varargin{1};
else
    r = 1;
end
A = deg2rad(A);
B = deg2rad(B);
dLat = A(:,1) - B(:,1);
dLong = A(:,2) - B(:,2);
res = 2 * r * asin(sqrt((sin(dLat/2)).^2 ...
    + cos(A(:,1)) .* cos(B(:,1)) .* (sin(dLong/2)).^2));

end