function [P,Pind] = perpLine(Q,Line)
% Finds the most perpendicular line PQ passing the point Q
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
if length(indmin)>1
    disp('Warning: the index is not unique')
    P = Line(:,indmin(1));
    Pind = indmin(1);    
    % Check if the given point was on the line
    xInd = find(Line(1,:)==Q(1)); % finds if x coord matches
    if Line(2,xInd)==Q(2) % if the corresponding y coord matches
        P = Q;
        Pind = xInd; % decide the query point is the "most perpendicular"
    end
else
    P = Line(:,indmin(1));
    Pind = indmin(1);
end