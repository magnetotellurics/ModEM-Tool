function [E,T,Grid,Modes] = rdExyz(cfile, fNum, mNum)
%  Usage [E,T,Grid,Modes] = rdExyz(cfile, fNum, mNum);
%  Reads in Ex, Ey, Ez; Periods, Grid, and list
%   of mode names from a standard fortran binary 
%  electric field for a given frequency and mode number
%  file cfile and fNum = frequency number; mNum = mode number
% NB Fortran stable- updated to always output integer*4  - so if this
% doesn't work for you anymore, fix recLengthForm [A. Kelbert, July 2020]

% first check the byte length of the record header for a fortran sequential
%     binary file
byteLength = fourOReight(cfile);
if byteLength == 4
    recLengthForm = 'integer*4';
elseif byteLength == 8
    recLengthForm = 'integer*8';
else
    error([cfile ' is not a Fortran sequential binary file'])
end
   
%  open file for reading
fid = fopen(cfile,'r');
% frewind(fid);
if(fid < 0) 
   error(['error: cannot open file ',cfile]);
end

% Generic Header
% First record: nPer, nMode, Nx, Ny, Nz, NzAir
% the first record is the generic header and all
% are integers
% skip past record length for record #1 
fseek(fid,byteLength,'cof');
%  read in version character string
cVersion = fread(fid,20,'char');
%  read in nPer, nMode, Nx, Ny, Nz, NzAir
N = fread(fid,6,'integer*4');
nPer = N(1);
nMode = N(2);
Nx = N(3);Ny = N(4); Nz = N(5);
NzAir = N(6); 
origin = fread(fid,3,'float64');
Rot = fread(fid,1,'float64');
%  skip past record length for record #1 
fseek(fid,2*byteLength,'cof');
dx = fread(fid,Nx,'float64');
fseek(fid,2*byteLength,'cof');
dy = fread(fid,Ny,'float64');
fseek(fid,2*byteLength,'cof');
dz = fread(fid,Nz,'float64');
fseek(fid,byteLength,'cof');
Grid = struct('dx',dx,'dy',dy,'dz',dz,'NzAir',NzAir, ...
    'origin', origin, 'rotation', Rot);

% we have skipped the main headers, including dx, dy, dz 
mainHeader = ftell(fid);

% first, determine the size of data sets we will
% be dealing with and then we decide how much data
% needs to be skipped
% skipping through the frequency header (36 + 8 = 44 bytes)
lrec = fread(fid,1,recLengthForm);
fseek(fid,byteLength+lrec,'cof');

% determining the size of Ex
lrec = fread(fid,1,recLengthForm);
% skipping through Ex data
fseek(fid,byteLength+lrec,'cof');

% determining the size of Ey
lrec = fread(fid,1,recLengthForm);
% skipping through Ey data
fseek(fid,byteLength+lrec,'cof');

% determining the size of Ez
lrec = fread(fid,1,recLengthForm);
fseek(fid,byteLength+lrec,'cof');
offset = ftell(fid);
% we have determined the length of Ex, Ey, and Ez

%   get list of periods
T = zeros(1,nPer);
for k = 1:nPer
   nSkip = (k-1)*nMode*(offset-mainHeader)+mainHeader+byteLength;
   fseek(fid,nSkip,'bof');
   omega = fread(fid,1,'float64');
   T(k) = 2*pi/omega;
end

for k = 1:nMode
   nSkip = (k-1)*(offset-mainHeader)+mainHeader+16+byteLength;
   fseek(fid,nSkip,'bof');
   Modes{k} = char(fread(fid,20,'char')');
end

%   get list of modes

whereTo = (fNum-1)*nMode+(mNum-1);
% main header of 52 bytes + Grid overhead is only skipped once
% offset-mainHeader is the frequency header and the 
% corresponding data
skipData = whereTo*(offset-mainHeader);
skipData = skipData + mainHeader;
% rewind and reposition the position indicator for the file
% frewind(fid);
fseek(fid, skipData, 'bof');

% Specific Header
% Second record: header for that frequency and 
% mode (36 + 8 = 44 bytes)
% skip past record length for record #2
lrec = fread(fid,1,recLengthForm);
% read the angular frequency, omega
inOmega = fread(fid,1,'real*8');
% read the frequency number
ifreq = fread(fid,1,'integer*4');
% read the mode number
imode = fread(fid,1,'integer*4');
% read whichMode
whichMode = char(fread(fid,20,'char'));
% skip past record length for record #2
lrec = fread(fid,1,recLengthForm);

% check whether we have the correct frequency and
% the mode
if (ifreq ~= fNum)
    fprintf(1,'Mismatch in Frequency\n')
end
if (imode ~= mNum)
    fprintf(1,'Mismatch in Mode\n')
end
%  Thrid record: Ex, Ey, and Ez
%  skip past record length for record #3 
lrec = fread(fid,1,recLengthForm);
%  read in Ex (NOTE: no complex IO in matlab)
NT = Nx*(Ny+1)*(Nz+1);
TEMP = fread(fid,[2,NT],'real*8');
Ex = reshape(complex(TEMP(1,:),TEMP(2,:)),[Nx,Ny+1,Nz+1]);
fseek(fid,2*byteLength,'cof');
NT = (Nx+1)*Ny*(Nz+1);
TEMP = fread(fid,[2,NT],'real*8');
Ey = reshape(complex(TEMP(1,:),TEMP(2,:)),[Nx+1,Ny,Nz+1]);
fseek(fid,2*byteLength,'cof');
NT = (Nx+1)*(Ny+1)*Nz;
TEMP = fread(fid,[2,NT],'real*8');
Ez = reshape(complex(TEMP(1,:),TEMP(2,:)),[Nx+1,Ny+1,Nz]);


E = struct('x',Ex,'y',Ey,'z',Ez,'mode',whichMode,'omega',inOmega);
