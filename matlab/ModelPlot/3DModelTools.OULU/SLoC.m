function res = SLoC(A,B,C)
% The Spherical Law of Cosines
%   cos(alpha) = [cos(c) - cos(a)cos(b)] / sin(a)sin(b)
%   where
%       alpha = the spherical angle between AB and AC
%       a = the Great-circle distance of AB
%       b = the Great-circle distance of AC
%       c = the Great-circle distance of BC

a = GCD(A,B); b =  GCD(A,C); c = GCD(B,C);
res = (cos(c) - cos(a).*cos(b)) ./ (sin(a).*sin(b));

end

