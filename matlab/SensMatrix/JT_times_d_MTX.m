function [g] = JT_times_d_MTX(J,d);
%  variant on JT_times_d ... does not sum over
%   frequencies, and also returns "imaginary part",
%   for transmitters in input data vector d which are complex
%   output is a cell array of size (3,nTx);
%   g{3,:) is character string "Cmplx" or "Real";
%   g{1:2,:} are model parameter structures; only g{1,k} is
%   used if g{3,k} == 'Real'

nTx = length(d);
nJ = length(J{1});
if nJ ~= nTx
   fprintf(1,'%s\n',...
        'Error: data vector and sensitivity size not consistent')
   return
end

g = cell(3,nTx);
kk = 0;
for k = 1:nTx
   g{1,k} = J{2}{1}; g{1,k}.v = 0;
   nSite = length(d{k}.siteLoc);
   if J{1}{k}.Cmplx
      g{2,k} = g{1,k}; 
      g{3,k} = 'Cmplx'; 
      for j = 1:nSite
         kk = kk + 1;
         g{1,k}.v = g{1,k}.v + real(d{k}.Z(j))*J{2}{kk}.v;
         g{2,k}.v = g{2,k}.v + imag(d{k}.Z(j))*J{2}{kk}.v;
         kk = kk + 1;
         g{1,k}.v = g{1,k}.v + imag(d{k}.Z(j))*J{2}{kk}.v;
         g{2,k}.v = g{2,k}.v - real(d{k}.Z(j))*J{2}{kk}.v;
      end
   else
      g{3,k} = 'Real'; 
      for j = 1:nSite
         kk = kk + 1;
         g{1,k}.v = g{1,k}.v + d{k}.Z(j)*J{2}{kk}.v;
      end
   end
end
