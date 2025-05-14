function array = objArray(varargin)
% creates on array of objects

L = length(varargin);
for i = 1:L
    array(i) = varargin{i};
end