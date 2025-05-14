function [P,Pind] = dotProdMin(Q,Line)
% Finds the most perpendicular PQ line passing the point Q
% and one of the grid points P at the Line.
%
% INPUT
% Q     - 2d point
% Line  - 2xN vector of a line
%
% OUTPUT
% P     - minimiser
% Pind  - index of minimiser

% the number of grid points
N = length(Line(1,:));

% perpendicular line has inverted slope with minus sign
slope = (Line(2,1)-Line(2,N))/(Line(1,1)-Line(1,N));
invSlope = 1/slope;
slopes = zeros(1,N);

for i = 1:N
    slopes(i) = (Q(2)-Line(2,i))/(Q(1)-Line(1,i));
end
slopeDiff = abs(slopes+invSlope);
indmin = find(slopeDiff==min(slopeDiff));
P = Line(:,indmin);
Pind = indmin;