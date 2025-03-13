function sOut = CovMult(sIn,smthParams,nTimes)
%  Usage  : sout = CovMult(sIn,smthParams)
%  Input : sIn = input conductivity structure
%          smthParams = structure with paramters to control smoothing
%          nTimes = optional argument: number of times to apply the smoother
%		(default is 1)
%  Output : sOut = smoothed output conductivity structure

if nargin < 3
   nTimes = 1;
end

grid = smthParams.grid;
%   construct horizontal smoothing weights
Dy = [grid.Dy(1) grid.Dy grid.Dy(end)];
Cy = (Dy(1:end-1)+Dy(2:end))/2;
Cmin = min(Cy);
alphaMax = Cmin^2/6;
alphaMin = smthParams.rho*alphaMax;
z = [0 cumsum(grid.Dz(grid.Nza+1:end))];
z = (z(1:end-1)+z(2:end))/2;
zWt = min(1,z/smthParams.zMax);
alpha = zWt*alphaMax+(1-zWt)*alphaMin;
%  set up weights for smoothing horizontally ... initially based on diffusion
Wm = (2./(Cy(1:end-1).*(Cy(1:end-1)+Cy(2:end))))'*alpha;
Wm(1,:) = 0;
Wp = (2./(Cy(2:end).*(Cy(1:end-1)+Cy(2:end))))'*alpha;
Wp(end,:) = 0;
W0 = 1 - Wm-Wp;
%   but let's symmetrize!
Wm = (Wm(2:end,:)+Wp(1:end-1,:))/2;
Wp = Wm;
%  set up weights for smoothing vertically
Zm = .25*ones(size(sIn.v));
Zp = Zm;
Zm(:,1) = 0;
Zp(:,end) = 0;
Z0 = 1-Zm-Zp;
Zm = Zm(:,2:end);
Zp = Zp(:,1:end-1);
%initialize output
sOut = sIn;
nSmooth = smthParams.nOuter*nTimes;
%  outer loop: smooth horizontally, then vertical nSmooth times
for k = 1:nSmooth
    for j = 1:smthParams.nYinner
        sOut.v = W0.*sIn.v;
        sOut.v(2:end,:) = sOut.v(2:end,:)+Wm.*sIn.v(1:end-1,:);
        sOut.v(1:end-1,:) = sOut.v(1:end-1,:)+Wp.*sIn.v(2:end,:);
        sIn = sOut;
    end
    % skip vertical smoothing on last step to insure symmetry
    if k < nSmooth
        sIn = sOut;
        for j = 1:smthParams.nZinner
            sOut.v = Z0.*sIn.v;
            sOut.v(:,2:end) = sOut.v(:,2:end)+Zm.*sIn.v(:,1:end-1);
            sOut.v(:,1:end-1) = sOut.v(:,1:end-1)+Zp.*sIn.v(:,2:end);
        end
        sIn = sOut;
    end
end
