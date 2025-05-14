function ind = findNearest(Q,grid)
% Finds the index of the nearest value to Q in the given grid
%
% INPUT
% Q vector
% grid vector
%
% OUTPUT
% ind: nearest index

NQ = length(Q);
ind = zeros(NQ,1);

for i = 1:NQ;
    % calculate d^2 to each point in a grid
    delta = grid-Q(i);
    d2 = delta.^2;
    
    % find the least d^2 value and it's index
    indmin = find(d2==min(d2));
    % save the index, for the case of multiple index
    ind = indmin;
    %ind(i) = indmin(1);
%     if length(indmin) > 1
%         disp('Warning: the nearest index was not unique!');
%     end
end