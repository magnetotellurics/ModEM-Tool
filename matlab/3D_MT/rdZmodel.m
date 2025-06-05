%************************************************************************Z
% readZ(cfile, fNum): Reads impedances from all the sites from a 3D forward modeling
% USAGE:  [Z,loc,T,nPer] = rdZmodel(cfile, fNum);
% cfile = input file name
% fNum = frequency number (input)
%  Z(2,2,nSites) = array of impedances 
%        (plus vertical field TFs if in file ... then Z(3,2,nSites))
%  loc(3,nSites) = array of locations

function [Z,loc,T,nPer] = readZ(cfile, fNum)

%  open file for reading
fid = fopen(cfile,'r');

% Generic Header for impedance file
% First record: nPer, nSites, ifBzTF
% the first record is the generic header and all
% are integers
% skip past record length for record #1 
lrec = fread(fid,1,'integer*4');
if lrec == 32
   %   this has version number in header
   %  read in version number
   fVer = fread(fid,20,'char');
end

%  read in nPer, nSites, ifBzTF
N = fread(fid,3,'integer*4');
nPer = N(1);
nSites = N(2);
ifBzTF = N(3);
if ifBzTF
   NTF = 3;
else
   NTF = 2;
end

%  skip past record length for record #1 
lrec = fread(fid,1,'integer*4');
% we have skipped the main header for an impedance 
% file (contains 12 + 8 = 20 bytes)
mainHeader = ftell(fid);

% Determine array sizes to decide how much data
% needs to be skipped
lrec = fread(fid,1,'integer*4');
status = fseek(fid,4+lrec,'cof');

% getting into the specifics of impedance data: iSite, xcCoord, ycCoord,
% zcCoord
% skipping through the specifics of the impedance data site (28 + 8 = 36
% bytes)
lrec = fread(fid,1,'integer*4');
status = fseek(fid,4+lrec,'cof');

lrec = fread(fid,1,'integer*4');
% skipping through Zxx (16 + 8 = 24 bytes)
status = fseek(fid,4+lrec,'cof');

lrec = fread(fid,1,'integer*4');
% skipping through Zxy (16 + 8 = 24 bytes)
status = fseek(fid,4+lrec,'cof');

lrec = fread(fid,1,'integer*4');
% skipping through Zyx (16 + 8 = 24 bytes)
status = fseek(fid,4+lrec,'cof');

lrec = fread(fid,1,'integer*4');
% skipping through Zyy (16 + 8 = 24 bytes)
status = fseek(fid,4+lrec,'cof');

if (ifBzTF)
   
    lrec = fread(fid,1,'integer*4');
    % skipping through Zzx (16 + 8 = 24 bytes)
    status = fseek(fid,4+lrec,'cof');

    lrec = fread(fid,1,'integer*4');
    % skipping through Zzy (16 + 8 = 24 bytes)
    status = fseek(fid,4+lrec,'cof');   
    
end    

offset = ftell(fid);
% we have determined the length of impedance data set for a given site and
% a given frequency

whereTo = (fNum-1)*nSites;
% main header of 12 + 8 = 20 bytes is only skipped once
% offset-mainHeader is the frequency header and the 
% corresponding impedance data
skipData = whereTo*(offset-mainHeader);
skipData = skipData + mainHeader;
% rewind and reposition the position indicator for the file
frewind(fid);
status = fseek(fid, skipData, 'bof');

k = 1;
%  Read Data
iSite = 0;
loc = zeros(3,nSites);
Z = zeros(NTF,2,nSites);

while (k <= nSites)
    
    % Specific Header
    % Second record: header for that frequency in impedance data (16 + 8 = 24 bytes) 
    % skip past record length for record #2
    lrec = fread(fid,1,'integer*4');
    % read the angular frequency, omega
    inOmega = fread(fid,1,'real*8');
    T = 2*pi/inOmega;
    % read the frequency number
    ifreq = fread(fid,1,'integer*4');
    % skip past record length for record #2
    lrec = fread(fid,1,'integer*4');

    % getting into the specifics of impedance data: iSite, xcCoord, ycCoord,
    % zcCoord (28 + 8 = 36 bytes)
    % skipping through the specifics of the impedance data (36 bytes)
    lrec = fread(fid,1,'integer*4');
    iSite = fread(fid,1,'integer*4');
    loc(:,iSite) = fread(fid,3,'real*8');
    lrec = fread(fid,1,'integer*4');

    lrec = fread(fid,1,'integer*4');
    % Zxx (16 bytes)
    temp = fread(fid,2,'real*8');
    Z(1,1,iSite) = complex(temp(1),temp(2));
    lrec = fread(fid,1,'integer*4');    

    lrec = fread(fid,1,'integer*4');
    % Zxy (16 bytes)
    temp = fread(fid,2,'real*8');
    Z(1,2,iSite) = complex(temp(1),temp(2));
    lrec = fread(fid,1,'integer*4'); 
    
    lrec = fread(fid,1,'integer*4');
    % Zyx (16 bytes)
    temp = fread(fid,2,'real*8');
    Z(2,1,iSite) = complex(temp(1),temp(2));
    lrec = fread(fid,1,'integer*4'); 
    
    lrec = fread(fid,1,'integer*4');
    % Zyy (16 bytes)
    temp = fread(fid,2,'real*8');
    Z(2,2,iSite) = complex(temp(1),temp(2));
    lrec = fread(fid,1,'integer*4'); 
    
    if (ifBzTF)
   
        lrec = fread(fid,1,'integer*4');
        % Zzx (16 bytes)
        temp = fread(fid,2,'real*8');
        Z(3,1,iSite) = complex(temp(1),temp(2));
        lrec = fread(fid,1,'integer*4'); 
        
        lrec = fread(fid,1,'integer*4');
        % Zzy (16 bytes)
        temp = fread(fid,2,'real*8');
        Z(3,2,iSite) = complex(temp(1),temp(2));
        lrec = fread(fid,1,'integer*4'); 
        k = k + 1;
        
    end    

end
