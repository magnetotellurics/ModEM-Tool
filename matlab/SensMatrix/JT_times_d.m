function [mOut] = JT_times_d(J,d)
%  Usage : [mOut] = JT_times_d(J,d);
% multiply by transpose of sensitivity matrix (stored as 
%     cell array of data/model space structures)

nTx = length(d);
nJ = length(J{1});
if nJ ~= nTx  | lengthDat(J{1}) ~= lengthDat(d)
   fprintf(1,'%s\n',...
	'Error: data vector and sensitivity size not consistent')
   return
else
end
kk = 0;
mOut = J{2}{1};
mOut.v = zeros(size(mOut.v));
for k = 1:nTx
    nSite = length(d{k}.siteLoc);
    for j = 1:nSite
       kk = kk + 1;
       mOut.v = mOut.v + real(d{k}.Z(j))*J{2}{kk}.v;
       if J{1}{k}.Cmplx
          kk = kk + 1;
          mOut.v = mOut.v + imag(d{k}.Z(j))*J{2}{kk}.v;
       end
    end
end
