function [I,J,X,Y] = sortLoc(loc);
%   given a list of locations on a regular grid
%   construct (i,j) index arrays for lines of constant x,y
%  Usage: [I,J,X,Y] = sortLoc(loc);

X = loc(1,1); Y = loc(2,1);
Nloc = length(loc);
for k = 2:Nloc
   if(~any(loc(1,k) == X))
      X = [X loc(1,k)];
   end
   if(~any(loc(2,k) == Y))
      Y = [Y loc(2,k)];
   end
end
I = zeros(Nloc,1); J = zeros(Nloc,1);
for k = 1:length(X)
   I(find(loc(1,:)==X(k))) = k;
end
for k = 1:length(Y)
   J(find(loc(2,:)==Y(k))) = k;
end
