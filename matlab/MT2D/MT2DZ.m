classdef MT2DZ < DataVector
    % base class for storing 2D MT data; follows structure used in 
    % original non-object-oriented "dataSpace" routines; idea is to
    % quickly develop a data-space object that can later be replaced by
    % something more logical

   properties
   %   NTX  %    Integer giving number of transmitters
   %   Nd   %    Integer: total number of real data
   %   d    %    cell array, one cell for each transmitter
   %   normalized = false   %   logical variable: is data vector normalized
         %     by standard deviations
   end

   methods

       %*******************************************************************
       function obj = MT2DZ(ntx)
       % class constructor 
          if nargin > 0
             obj.NTX  = ntx;
             obj.d = cell(ntx,1);
          else
             obj.NTX = 0;
          end
       end
       %*******************************************************************
       function obj = regTemplate(obj,m,OPTIONS)
       %   fill  in data object template  ... 
       %                           ... First site locations
           if nargin==2 
               
               OPTIONS = struct('CellCenters',true,'CellLim',8,...
                   'Nsites',0,'SiteLims',[0,0],...
                   'Nperiods',20,'PeriodRange',[10,10000],...
                   'Seafloor',true,'MODE','JT');
           end
           y = cumsum([0 m.grid.Dy]);
           yCtr = (y(1:end-1)+y(2:end))/2;
           z = cumsum([0 m.grid.Dz]);
           if OPTIONS.CellCenters 
               sitesY = yCtr(OPTIONS.CellLim:end-OPTIONS.CellLim-1);
               sitesZ = z(m.grid.Nza+1)*ones(size(sitesY));
               Nsites = length(sitesY);
           else
               error('Not coded for anything but cell centers')
           end
           
           %   ocean floor sites
           if OPTIONS.Seafloor
               for k = 1:Nsites
                   j =  find(sitesY(k) < y,1);
                   %    assume any surface cell exceeding 3 S/m is ocean
                   if(m.v(j,1)> log(3))
                       %   look for bottom of column of ocean cells
                       kk  =  find(m.v(j,:)> log(3),1,'last');
                       sitesZ(k) = z(m.grid.Nza+kk+1);
                   end
               end
           end
           siteLoc = [sitesY' sitesZ'];
           siteChar  =  [];
           for k = 1:Nsites
               siteChar = strvcat(siteChar,num2str(k));
           end
           Z = zeros(Nsites,1);
           Zerr = ones(Nsites,1);
           
           % Next periods
           T2 = OPTIONS.PeriodRange(2);
           T1 = OPTIONS.PeriodRange(1);
           nPer = OPTIONS.Nperiods;
           dt = (log10(T2/T1))/(nPer-1);
           periods = 10.^(log10(T1):dt:log10(T2));

           if strcmp(OPTIONS.MODE,'JT')
               %  joint TE/TM (at all sites)
               %   first all TE
               obj.NTX = 2*nPer;
               obj.d = cell(2*nPer,1);
                for j=1:nPer
                    obj.d{j} = struct('T',periods(j),...
                       'Cmplx',1,...
	               'Mode','TE',...
	               'siteLoc',siteLoc,...
                   'siteChar',siteChar,...
	               'Z',Z,'Zerr',Zerr);
               end 
               %   then all TM
               for j=nPer+1:2*nPer
                    obj.d{j} = struct('T',periods(j-nPer),...
                       'Cmplx',1,...
	               'Mode','TM',...
	               'siteLoc',siteLoc,...
                   'siteChar',siteChar,...
	               'Z',Z,'Zerr',Zerr);
               end
           else
           % only one of TE and TM ...
               obj.NTX = nPer;
               obj.d = cell(nPer,1);
               for j=1:nPer
                   obj.d{j} = struct('T',periods(j),...
                       'Cmplx',1,...
	               'Mode',OPTIONS.MODE,...
	               'siteLoc',siteLoc,...
                   'siteChar',siteChar,...
	               'Z',Z,'Zerr',Zerr);
               end
           end
           obj.Nd = obj.NTX*Nsites;
       end
       %**********************************************************************
       function [obj] = readVec(obj,cfile)
       %
       %  Usage: [obj,status] = readVec(obj,cfile)
       %
       %   reads data vector from impedance file, puts contents into
       %     empty (but already created) MT2DZ (2D MT data vector) object
           fid = fopen(cfile,'r');

           %  first read the header
           record = textscan(fid,'%*[^:]%*c%[^\n]',3);
           %info = record{1}{1};
           %units = record{1}{2};
           %isign = str2num(record{1}{3});

           %  record # 1: number of transmitters
           nTx = fscanf(fid,'%d\n',1);
           obj.NTX =  nTx;
           obj.Nd = 0;

           for j = 1:nTx
               %   for each period/mode ....

               %  record # 1 in transmitter block: period, mode, number of sites
               record = textscan(fid,'%f %2c %d\n',1);
               T = record{1};
               MODE = record{2};
               nSites = double(record{3});
               nComp = 2;
               obj.Nd = obj.Nd+nComp*nSites;

               %  record # 2 in transmitter block: site locations
               record = fscanf(fid,'%f',[nSites,2]);
               siteLoc = record;
               record = fscanf(fid,'\n');

               %  read and ignore the header of the data block
               fgetl(fid);

               %  records # 3 & #4 in transmitter block: data and error bars
               siteChar = '';
               Z = zeros(nSites,floor(nComp/2));
               Z = Z+1i*Z;
               Zerr = Z;
               for k = 1:nSites
                   idx = 1:2:nComp-1;
                   record = fscanf(fid,'%s',1);
                   siteChar = strvcat(siteChar, record);
                   record = fscanf(fid,'%f',nComp);
                   Z(k,:) = record(idx)+1i*record(idx+1);
                   record = fscanf(fid,'%f',nComp);
                   Zerr(k,:) = record(idx);
               end

               obj.d{j} = struct('T',T,'Cmplx',1,...
                   'Mode',MODE,'siteLoc',siteLoc,'siteChar',siteChar,'Z',Z,'Zerr',Zerr);

               clear siteChar
           end
           fclose(fid);
       end
       %*******************************************************************
       function obj = SetDataTX(obj,itx,T,MODE,siteLoc,Z,Zerr)
       %   add data structure for one transmitter

          Cmplx = 1-all(isreal(Z));

          [n1,~] = size(siteLoc);
          if n1 == 2
              siteLoc = siteLoc';
          end
          [n1,~] = size(Z);
          if n1 == 1
              Z = Z.';
          end
          [n1,~] = size(Zerr);
          if n1 == 1
              Zerr = Zerr';
          end
          obj.d{itx} = struct('T',T,'Cmplx',Cmplx,...
           'Mode',MODE,'siteLoc',siteLoc,'Z',Z,'Zerr',Zerr);
          if itx > obj.NTX
              obj.NTX = itx;
          end
       end
       %*******************************************************************
       function [status] = writeVec(obj,cfile,info,units,isign)
       %  Usage:  [status] = writeZ_2D(cfile,allData,info,units,isign);
       %   arguments info, units and isign are optional.
       %   write contents of cell array allData to file
       %   cfile.  There is one cell per period; each
       %   cell contains all information necessary to define
       %   data (locations, values, error standard dev) for
       %   each period (transmitter)
       %   NB  Assuming complex data
           
           fid = fopen(cfile,'w');
           
           %  description <= 80 char in length
           if nargin < 3
               info = 'Synthetic 2D MT data written in Matlab';
           end
           
           %  assume SI units by default; alternative is [mV/km]/[nT]
           if nargin < 4
               units = '[V/m]/[T]';
           end
           
           %  sign convention is -1 by default
           if nargin < 5
               isign = -1;
           end
           
           %  header: three lines followed by an empty line
           fprintf(fid,'Description: %s\n',info);
           fprintf(fid,'Units: %s\n',units);
           fprintf(fid,'Sign convention: %d\n\n',isign);
           
           %  record # 1: number of transmitters
           status = fprintf(fid,'%d\n',obj.NTX);
           
           for j = 1:obj.NTX
               
               %   for each period/mode ....
               
               %  record # 1 in transmitter block: period, mode, number of sites
               T = obj.d{j}.T;
               MODE = obj.d{j}.Mode;
               nSites = size(obj.d{j}.siteLoc);
               nComp = 2;
               fprintf(fid,'%12.6E %5s %5d\n',T,MODE,nSites(1));
               
               %  record # 2 in transmitter block: site locations
               siteLoc = obj.d{j}.siteLoc;
               for i = 1:2
                   for k = 1:nSites
                       fprintf(fid,'%12.3f',siteLoc(k,i));
                   end
                   fprintf(fid,'\n');
               end
               
               %  header of the data & error block: will in future be more informative
               status = fprintf(fid,'%12s Re %12s Im\n',' ',' ');
               
               %  records # 3 & #4 in transmitter block: data and error bars
               for k = 1:nSites
                   if isfield(obj.d{j},'siteChar')
                       siteChar = obj.d{j}.siteChar(k,:);
                   else
                       siteChar = num2str(k);
                   end
                   fprintf(fid,'%10s',siteChar);
                   Zri = obj.d{j}.Z;
                   if obj.d{j}.Cmplx
                       for i = 1:nComp/2
                           fprintf(fid,'%15.6E %15.6E',real(Zri(k,i)),imag(Zri(k,i)));
                       end
                   else
                       for i = 1:nComp
                           fprintf(fid,'%15.6E',Zri(k,i));
                       end
                   end
                   fprintf(fid,'\n%10s',' ');
                   Zerr = obj.d{j}.Zerr;
                   if obj.d{j}.Cmplx
                       for i = 1:nComp/2
                           fprintf(fid,'%15.6E %15.6E',Zerr(k,i),Zerr(k,i));
                       end
                   else
                       for i = 1:nComp
                           fprintf(fid,'%15.6E',Zerr(k,i));
                       end
                   end
                   fprintf(fid,'\n');
               end
           end
       end
       %**************************************************************************
       function [dOut,dErr,normalized] = ExtractVec(obj)
       %  Makes a standard real vector out of an MT2DZ object
       %  Also returns data standard error as a vector
       %
       %  Usage: [dOut,dErr] = ExtractVec(obj);
       %         [dOut,dErr,normalized] = ExtractVec(obj);
       %    optional third argument returns the normalized attribute stored
       %    in the input data vector

          dOut = [];
          for k = 1:obj.NTX
             if obj.d{k}.Cmplx
                Zr = real(obj.d{k}.Z);
                Zi = imag(obj.d{k}.Z);
                Z = zeros(2*length(Zr),1);
                Z(1:2:end) = Zr;
                Z(2:2:end) = Zi;
             else
                Z = real(obj.d{k}.Z);
             end
             dOut = [dOut; Z];
          end
          if nargout >=2
             dErr = [];
             for k = 1:obj.NTX
                nsta = length(obj.d{k}.Zerr);
                if obj.d{k}.Cmplx
                   Zerr = zeros(1,2*nsta);
                   Zerr(1:2:end-1) = obj.d{k}.Zerr;
                   Zerr(2:2:end) = obj.d{k}.Zerr;
                else
                   Zerr = obj.d{k}.Zerr;
                end
                dErr = [dErr; Zerr'];
             end
          end
          if nargout == 3
              normalized = obj.normalized;
          end
       end
       %**********************************************************************
       function objOut = SetVec(objIn,dIn,normalized)
       %  Usage: objOut = SetVec(objIn,dIn);
       %  copies the contents of an simple vector object into an existing MT2DZ
       %    object, keeping sizes, order, station coordinates the same (size of
       %    input vector and MT2DZ object must be compatible)
       %   optional argument "normalized" is used to set attribute in
       %   output vector ... does not actually result in any 
       %   normalization/unnormalization inside routine
       
          objOut = objIn;
       
          if nargin == 3
              objOut.normalized = normalized;
          end

          ii = 1;
          for k = 1:objIn.NTX
             nSites = length(objIn.d{k}.siteLoc);
             for j = 1:nSites
                if objOut.d{k}.Cmplx
                   objOut.d{k}.Z(j) = dIn(ii)+1i*dIn(ii+1);
                   ii = ii + 2;
                else
                   objOut.d{k}.Z(j) = dIn(ii);
                   ii = ii + 1;
                end
             end
          end
       end
       %*******************************************************************
       function  obj = zeroVec(obj)  
       %  zeros the data values in an existing data vector object
           for k = 1:obj.NTX
               obj.d{k}.Z = zeros(size(obj.d{k}.Z));
           end   
       end    
       %*******************************************************************
       function obj = plus(obj1,obj2)     
       %  Implements addition for data vector objects ... simple at present,
       %    just copies metadata (and error bars) from first object
             %   so far error bars  are just copied from d1;
       %
       %  Usage : d = plus(d1,d2);
       %          d = d1+d2;

          if obj1.NTX ~= obj2.NTX
             error('Error: data vector objects not compatible')
          end
          
          if obj1.normalized ~= obj2.normalized
              error('Error: adding normalized and unnormalized data vectors')
          end

          obj = obj1;
          for k = 1:obj1.NTX
             obj.d{k}.Cmplx = obj1.d{k}.Cmplx | obj2.d{k}.Cmplx;
             obj.d{k}.Z = obj1.d{k}.Z + obj2.d{k}.Z;
          end
       end
       %*******************************************************************
       function obj = add1Tx(obj,j,c,d1)     
       %  adds c*d1, partial data vector for a single transmitter to data
       %     vector object
       %  Usage : obj = add1Tx(obj,k,c,d1) ;
          obj.d{j}.Z = obj.d{j}.Z + c*d1.d{j}.Z;
       end
       %*******************************************************************
       function obj = mult1Tx(obj,j,c)     
       %  multiplies data vector for a single transmitter (j) by a constant
       %  c
       %  Usage : obj = mult1Tx(obj,k,c) ;
          obj.d{j}.Z = c*obj.d{j}.Z;
       end
       %*******************************************************************
       function obj = minus(obj1,obj2)     
       %  Implements subtraction for data vector objects ... simple at present,
       %    just copies metadata (and error bars) from first object
       %   so far error bars  are just copied from d1;
       %
       %  Usage : d = plus(d1,d2);
       %          d = d1+d2;

          if obj1.NTX ~= obj2.NTX
             error('Error: data vector objects not compatible')
          end
       
          if obj1.normalized ~= obj2.normalized
              error('Error: subtracting normalized and unnormalized data vectors')
          end

          obj = obj1;
          for k = 1:obj1.NTX
             obj.d{k}.Cmplx = obj1.d{k}.Cmplx | obj2.d{k}.Cmplx;
             obj.d{k}.Z = obj1.d{k}.Z - obj2.d{k}.Z;
          end 
       end
       %*******************************************************************
       function obj = mtimes(c,obj1)     
       %  Implements scalar multiplication for data vector objects ... simple at present,
       %    c must be double, must be on the left in the expression obj = c*obj1
       %
       %  Usage : d = c*d1
       %          d = mtimes(c,d1);

          if ~isa(c,'double')
             error('Error: can only multiply data vectors by doubles')
          end

          obj = obj1;

          for k = 1:obj1.NTX
             obj.d{k}.Cmplx = obj1.d{k}.Cmplx | ~isreal(c);
             obj.d{k}.Z = c*obj1.d{k}.Z;
             %  make any changes to error bars explicit!
             %obj.d{k}.Zerr = abs(c)*obj1.d{k}.Zerr;
          end 
       end
       %*******************************************************************
       function obj = linComb(c1,obj1,c2,obj2)
       %  Usage : [d] = linCombDat(d1,c1,d2,c2);
       %   computes d = c1*d1+c2*d2 for two data vector object
       %    so far error bars  are just copied from d1;

          if obj1.NTX ~= obj2.NTX
             error('Error: data vector objects not compatible')
          end
          if obj1.normalized ~= obj2.normalized
              error('Error: linear  combination of normalized and unnormalized data vectors')
          end
          obj = obj1;
          for k = 1:obj.NTX
             obj.d{k}.Cmplx = obj1.d{k}.Cmplx | obj2.d{k}.Cmplx;
             obj.d{k}.Z = c1*obj1.d{k}.Z + c2*obj2.d{k}.Z;
          end
       end
       %*******************************************************************
       function n = lengthDat(obj)
       %  Total length of a data vector object
       %
       %  Usage [n] = lengthDat(d);
          n = 0;
          for k = 1:obj.NTX
             if obj.d{k}.Cmplx
                n = n + 2*length(obj.d{k}.siteLoc);
             else
                n = n + length(obj.d{k}.siteLoc);
             end
          end
       end
       %******************************************************************
       function  obj = Normalize(obj,Ntimes)
       %  Normalizes a data vector by dividing by error standard deviation
       %      Ntimes is 1 to divide by standard deviation
       %                2 to divide by variance
       %             This is optional, default is 1
       %   NO MODIFICATION TO ERROR STANDARD DEVIATION STORED IN OBJECT

          if nargin == 1
             Ntimes = 1;
          end
          if obj.normalized
              fprintf(1,'%s\n','Warning: data vector already normalized')
          else
             for k = 1:obj.NTX
                if Ntimes == 1
                   obj.d{k}.Z = obj.d{k}.Z./obj.d{k}.Zerr;
                else
                   obj.d{k}.Z = obj.d{k}.Z./(obj.d{k}.Zerr).^2;
                end
             end
             obj.normalized = true;
          end
       end
       %******************************************************************
       function  obj = UnNormalize(obj,Ntimes)
       %  "UnNormalizes" a data vector by multiplying by error standard deviation
       %      Ntimes is 1 to multiply by standard deviation
       %                2 to multiply by variance
       %             This is optional, default is 1
       %   NO MODIFICATION TO ERROR STANDARD DEVIATION STORED IN OBJECT
       
          if nargin == 1
             Ntimes = 1;
          end
          if ~obj.normalized
              fprintf(1,'%s\n','Warning: data vector already UnNormalized')
          else
             for k = 1:obj.NTX
                if Ntimes == 1
                   obj.d{k}.Z = obj.d{k}.Z.*obj.d{k}.Zerr;
                else
                   obj.d{k}.Z = obj.d{k}.Z.*(obj.d{k}.Zerr).^2;
                end
             end
             obj.normalized = false;
          end
       end
       %******************************************************************
       function CdInv = InvDataErr(obj)
       % Makes a standard real vector containing the inverse of the data
       % error standard deviation from an input impedance data vector object
       %
       % Usage:  CdInv = InvDataErr(dIn)
           CdInv = [];
           for k = 1:obj.NTX
               Zr = real(obj.d{k}.Zerr);
               if obj.d{k}.Cmplx
                   % errors for real and imaginary part are assumed equal
                   Cd = zeros(length(Zr)*2,1);
                   Cd(1:2:end) = Zr;
                   Cd(2:2:end) = Zr;
               else
                   Cd = Zr;
               end
               CdInv = [CdInv; Cd];
           end
           CdInv = 1./CdInv;
       end
       %*******************************************************************
       function  obj = MultCdInv(obj,objErr)
       %  Normalizes a data vector by dividing by error standard deviation
       %   whether the object is already normalized or not
       %   if a second argument is present errors from this data vector are
       %   used for the normalization, and these errors are copied into
       %    the output vector

           if nargin == 1
               for k = 1:obj.NTX
                   obj.d{k}.Z = obj.d{k}.Z./obj.d{k}.Zerr;
               end
           else
               for k = 1:obj.NTX
                   obj.d{k}.Z = obj.d{k}.Z./objErr.d{k}.Zerr;
                   obj.d{k}.Zerr = objErr.d{k}.Zerr;
               end
           end
           obj.normalized = true;
       end
       %*******************************************************************
       function Cd = ErrSD(dIn)
       % Makes a standard real vector containing the data error standard 
       % deviation extracted from input impedance data vector object
       %
       % Usage:  Cd = ErrSD(dIn);
           
           Cd = [];
           for k = 1:dIn.NTX
               Zr = real(dIn.d{k}.Zerr);
               if dIn.d{k}.Cmplx
                   % errors for real and imaginary part are assumed equal
                   Cd1 = zeros(length(Zr)*2,1);
                   Cd1(1:2:end) = Zr;
                   Cd1(2:2:end) = Zr;
               else
                   Cd1 = Zr;
               end
               Cd = [Cd; Cd1];
           end
       end
       %*******************************************************************
       function obj = addNoise(obj,relErr)
       %  Usage : addNoise(obj,relErr);
       %    adds gaussian random noise to data vector components
       %    with standard deviation relErr*abs(Z)
       %    Sets error standard deviation to appropriate value

          for k = 1:obj.NTX
              temp = obj.d{k}.Z;
              errScale = relErr*abs(temp);
              obj.d{k}.Z = obj.d{k}.Z + errScale.*(randn(size(obj.d{k}.Z))+...
                                         1i*randn(size(obj.d{k}.Z)));
              obj.d{k}.Zerr = errScale;
          end
       end
       %*******************************************************************
       function ip = dotCmplx(obj1,obj2)
       %   inner product of two (complex impedance)
       %    data objects, using error standard deviation 
       %    defined from d1 
       %  Usage :  ip = ipDat(d1,d2);

          if obj1.NTX ~= obj2.NTX
              error('incompatible data objects for inner product')
          end
          if obj1.normalized ~= obj2.normalized
              error('Error: inner product of normalized and unnormalized data vectors')
          end
          ip= 0.;
          if obj1.normalized
              for k = 1:obj1.NTX
                 ip = ip + obj1.d{k}.Z'*obj2.d{k}.Z;
              end 
          else
              for k = 1:obj1.NTX
                  ip = ip + (obj1.d{k}.Z./obj1.d{k}.Zerr)'* ...
                      (obj2.d{k}.Z./obj2.d{k}.Zerr);
              end
          end
       end
       %*******************************************************************
       function ip = dot(obj1,obj2)
       %   REAL inner product of two MT2DZ (possibly complex)
       %    data objects, using error standard deviation 
       %    defined from d1, if both objects are unnormalized 
       %  Usage :  ip = dot(d1,d2);

          if obj1.NTX ~= obj2.NTX
              error('incompatible data objects for inner product')
          end
          if obj1.normalized ~= obj2.normalized
              error('Error: inner product of normalized and unnormalized data vectors')
          end
          ip= 0.;
          if obj1.normalized
              for k = 1:obj1.NTX
                 ip = ip + obj1.d{k}.Z'*obj2.d{k}.Z;
              end 
          else
              for k = 1:obj1.NTX
                  ip = ip + (obj1.d{k}.Z./obj1.d{k}.Zerr)'* ...
                      (obj2.d{k}.Z./obj2.d{k}.Zerr);
              end
          end
          ip = real(ip);
       end
       %*******************************************************************
       function ip = dotMTX(obj1,obj2)
       %   REAL inner product of two (possibly complex impedance)
       %    data objects, using error standard deviation 
       %    defined from d1 (if both objects are unnormalized  ... 
       %    operates transmitter by transmitter, returns an array in ip
       %  Usage :  ip = dotMTX(d1,d2);

          if obj1.NTX ~= obj2.NTX
              error('incompatible data objects for inner product')
          end
          if obj1.normalized ~= obj2.normalized
              error('Error: inner product of normalized and unnormalized data vectors')
          end
          ip= zeros(obj1.NTX,1);
          if obj1.normalized
              for k = 1:obj1.NTX
                 ip(k) = obj1.d{k}.Z'*obj2.d{k}.Z;
              end 
          else
              for k = 1:obj1.NTX
                  ip(k) = (obj1.d{k}.Z./obj1.d{k}.Zerr)'* ...
                      (obj2.d{k}.Z./obj2.d{k}.Zerr);
              end
          end
          ip = real(ip);
       end
       %*******************************************************************
       function ip = dotCmplxNoCov(obj1,obj2)
       %   inner product of two (complex impedance)
       %    data objects, No  normalization by standard deviations
       %  Usage :  ip = ipDat(d1,d2);

          if obj1.NTX ~= obj2.NTX
              error('incompatible data objects for inner product')
          end
          ip= 0.;
          for k = 1:obj1.NTX
              ip = ip + (obj1.d{k}.Z)'*(obj2.d{k}.Z);
          end
       end
       %**************************************************************************
       function ip = dotNoCov(obj1,obj2)
       %   REAL inner product of two (possibly complex impedance)
       %    data objects, No normalization by standard deviations
       %  Usage :  ip = ipDat(d1,d2);
           
           if obj1.NTX ~= obj2.NTX
               error('incompatible data objects for inner product')
           end
           ip= 0.;
           for k = 1:obj1.NTX
               ip = ip + (obj1.d{k}.Z)'*(obj2.d{k}.Z);
           end
           ip = real(ip);
       end
       %**************************************************************************
       function objOut = UnitVec(objIn,MTX)
       %   normalizes data vector so that it has unit norm; if optional
       %   argument MTX is present and MTX = true each transmitter is
       %   normalized separately (to have unit norm).
       
           if nargin ==1
               MTX = false;
           end
           
           objOut = objIn;        
           if MTX
               for k = 1:objIn.NTX
                   objOut.d{k}.Z = objIn.d{k}.Z ...
                       /sqrt(objIn.d{k}.Z'*objIn.d{k}.Z);
               end        
           else
               norm = 0;
               for k = 1:objIn.NTX
                   norm = norm+objIn.d{k}.Z'*objIn.d{k}.Z;
               end  
               norm = sqrt(norm);
               for k = 1:objIn.NTX
                   objOut.d{k}.Z = objIn.d{k}.Z/norm;
               end        
           end
       end
       %***********************************************
       function diff  = relDiff(obj1,obj2)
           temp = obj1-obj2;
           diff = dot(temp,temp)/dot(obj1,obj1);
           diff = sqrt(diff);
       end
    end    %    methods
end     %   classdef
