classdef MT3DZ < DataVector
    % base class for storing 3D MT data; follows structure used in 
    % MT2DZ
    %
    %   Gary D. Egbert, 2010
    %   College of Oceanic and Atmospheric Sciences
    %   Oregon State University

   properties
   %   NTX  %    Integer giving number of transmitters
   %   Nd   %    Integer: total number of real data
   %   d    %    cell array, one cell for each transmitter
   %   normalized = 0   %   integer variable: normalization state of data
         %          vector
         %     by standard deviations
   end

   methods

       %*******************************************************************
       function obj = MT3DZ(ntx)
       % class constructor 
          if nargin > 0
             obj.NTX  = ntx;
             obj.d = cell(ntx,1);
          else
             obj.NTX = 0;
          end
       end
       %**********************************************************************
       function [obj] = readVec(obj,cfile)
       %
       %  Usage: [obj,status] = readVec(obj,cfile)
       %
       %   reads data vector from ASCII 3D impedance file, puts contents into
       %     empty (but already created) MT2DZ (2D MT data vector) object
       %     calls readZ_3D to fill in data structure  
       
           [obj.d] = readZ_3D(cfile,'[V/m]/[T]');
           obj.NTX = length(obj.d);
           obj.Nd = length(obj);
       end
       %*******************************************************************
       function [status] = writeVec(obj,cfile,info,units,isign)
       %  Usage:  [status] = writeVwc(obj,cfile,info,units,isign);
       %   arguments info, units and isign are optional.
       %   write 3D data vector object to file cfile, by calling writeZ_3D 
           
           %  description <= 80 char in length
           if nargin < 3
                info = 'matlab object oriented inversion';
           end
           
           %  assume SI units by default; alternative is [mV/km]/[nT]
           if nargin < 4
               units = '[V/m]/[T]';
           end
           
           %  sign convention is -1 by default
           if nargin < 5
               isign = -1;
           end
           [status] = writeZ_3D(cfile,obj.d,info,units,isign);
       end
       %**************************************************************************
       function [dOut,dErr,normalized] = ExtractVec(obj)
       %  Makes a standard real vector out of an MT3DZ object
       %  Also returns data standard error as a vector
       %
       %  Usage: [dOut,dErr] = ExtractVec(obj);
       %         [dOut,dErr,normalized] = ExtractVec(obj);
       %    optional third argument returns the normalized attribute stored
       %    in the input data vector

           dOut = zeros(obj.Nd,1);
           k1 = 1;
           for k = 1:obj.NTX
               [n1,n2] = size(obj.d{k}.Z);
               if obj.d{k}.Cmplx
                   k2 = k1+2*n1*n2-1;
                   Zr = real(obj.d{k}.Z);
                   Zi = imag(obj.d{k}.Z);
                   dOut(k1:2:k2-1) = reshape(Zr.',n1*n2,1);
                   dOut(k1+1:2:k2) = reshape(Zi.',n1*n2,1);
               else
                   k2 = k1+n1*n2-1;
                   dOut(k1:k2) = reshape(real(obj.d{k}.Z)',n1*n2,1);
               end
               k1 = k2+1;
           end

           if nargout >=2
               dErr = zeros(obj.Nd,1);
               k1 = 1;
               for k = 1:obj.NTX
                   [n1,n2] = size(obj.d{k}.Z);
                   if obj.d{k}.Cmplx
                       k2 = k1+2*n1*n2-1;
                       dErr(k1:2:k2-1) = reshape(obj.d{k}.Zerr',n1*n2,1);
                       dErr(k1+1:2:k2) = reshape(obj.d{k}.Zerr',n1*n2,1);   
                   else
                       k2 = k1+n1*n2-1;
                       dErr(k1:k2) = reshape(real(obj.d{k}.Zerr'),n1*n2,1);
                   end
                   k1 = k2+1;
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
           
           k1 = 1;
           for k = 1:objIn.NTX
               [n1,n2] = size(objIn.d{k}.Z);
               if objOut.d{k}.Cmplx
                   k2 = k1+2*n1*n2-1;
                   objOut.d{k}.Z = ...
                       reshape(dIn(k1:2:k2-1)+1i*dIn(k1+1:2:k2),n2,n1).';
               else
                   k2 = k1+n1*n2-1;
                   objOut.d{k}.Z = reshape(dIn(k1:k2-1),n2,n1).';
               end
               k1 = k2+1;
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
       function obj = add1Tx(obj,k,c,d1)     
       %  adds c*d1, partial data vector for a single transmitter to data
       %     vector object
       %  Usage : obj = add1Tx(obj,k,c,d1) ;
          obj.d{k}.Z = obj.d{k}.Z + c*d1.d{k}.Z;
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
       %   DON'T USE THIS -- USE length
       %
       %  Usage [n] = lengthDat(d);
          n = 0;
          for k = 1:obj.NTX
             [n1,n2] = size(obj.d{k}.Z); 
             if obj.d{k}.Cmplx
                n = n + 2*n1*n2;
             else
                n = n + n1*n2;
             end
          end
       end
       %*******************************************************************
       function n = length(obj)
       %  Total length of a data vector object
       %
       %  Usage [n] = length(d);
          n = 0;
          for k = 1:obj.NTX
             [n1,n2] = size(obj.d{k}.Z); 
             if obj.d{k}.Cmplx
                n = n + 2*n1*n2;
             else
                n = n + n1*n2;
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
             obj.normalized = Ntimes;
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
             Ntimes = obj.normalilzed;
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
             obj.normalized = Ntimes;
          end
       end
       %******************************************************************
       function CdInv = InvDataErr(obj)
       % Makes a standard real vector containing the inverse of the data
       % error standard deviation from an input impedance data vector object
       %
       % Usage:  CdInv = InvDataErr(dIn)
           CdInv = zeros(obj.Nd,1);
           k1 = 1;
           for k = 1:obj.NTX
               [n1,n2] = size(obj.d{k}.Z);
               if obj.d{k}.Cmplx
                   k2 = k1+2*n1*n2-1;
                   CdInv(k1:2:k2-1) = reshape(obj.d{k}.Zerr',n1*n2,1);
                   CdInv(k1+1:2:k2) = reshape(obj.d{k}.Zerr',n1*n2,1);
               else
                   k2 = k1+n1*n2-1;
                   CdInv(k1:k2) = reshape(obj.d{k}.Zerr',n1*n2,1);
               end
               k1 = k2+1;
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
       function obj = setErrScales(obj,relErr)
       %  Usage : seErrScales(obj,relErr);
       %    set error scales to relErr (fractional error; main diagonals 
       %    vertical field TFs treated differently.... this is QD
       
          for k = 1:obj.NTX
              temp = obj.d{k}.Z;
              errScale = relErr*abs(temp);
              [nSites,nComp] = size(temp);
              mainDiagScale = ones(nSites,1);
              for j = 1:nComp
                  if(strcmpi(obj.d{k}.compChar(j,:),'zxy') || ...
                      strcmpi(obj.d{k}.compChar(j,:),'zyx'));
                      mainDiagScale = mainDiagScale.*errScale(:,j);
                  end
              end
              mainDiagScale = sqrt(mainDiagScale);
              for j = 1:nComp
                  if (strcmpi(obj.d{k}.compChar(j,:),'zxx') || ...
                      strcmpi(obj.d{k}.compChar(j,:),'zyy'))
                      errScale(:,j) = mainDiagScale;
                  else
                      if (strcmpi(obj.d{k}.compChar(j,:),'tx ') || ...
                      strcmpi(obj.d{k}.compChar(j,:),'ty '))
                          errScale(:,j) = relErr*.25*ones(nSites,1);
                      end
                  end           
              end
              obj.d{k}.Zerr = errScale;
          end
       end
       %*******************************************************************
       function obj = addNoise(obj,relErr)
       %  Usage : addNoise(obj,relErr);
       %    adds gaussian random noise to data vector components
       %    with standard deviation relErr*abs(Z)
       %    Sets error standard deviation to appropriate value

          obj = setErrScales(obj,relErr);
          for k = 1:obj.NTX
              obj.d{k}.Z = obj.d{k}.Z + obj.d{k}.Zerr.*(randn(size(obj.d{k}.Z))+...
                                         1i*randn(size(obj.d{k}.Z)));
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
                 ip = ip + sum(sum(conj(obj1.d{k}.Z).*obj2.d{k}.Z));
              end 
          else
              for k = 1:obj1.NTX
                  ip = ip + sum(sum(conj(obj1.d{k}.Z./obj1.d{k}.Zerr).* ...
                      (obj2.d{k}.Z./obj2.d{k}.Zerr)));
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
                 ip = ip + sum(sum(conj(obj1.d{k}.Z).*obj2.d{k}.Z));
              end 
          else
              for k = 1:obj1.NTX
                  ip = ip + sum(sum(conj(obj1.d{k}.Z./obj1.d{k}.Zerr).* ...
                      (obj2.d{k}.Z./obj2.d{k}.Zerr)));
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
                 ip(k) = sum(sum(conj(obj1.d{k}.Z).*obj2.d{k}.Z));
              end 
          else
              for k = 1:obj1.NTX
                  ip(k) = sum(sum(conj(obj1.d{k}.Z./obj1.d{k}.Zerr).* ...
                      (obj2.d{k}.Z./obj2.d{k}.Zerr)));
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
              ip = ip + sum(sum(conj(obj1.d{k}.Z).*(obj2.d{k}.Z)));
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
               ip = ip + sum(sum(conj(obj1.d{k}.Z).*(obj2.d{k}.Z)));
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
                   temp = sum(sum(conj(objIn.d{k}.Z).*objIn.d{k}.Z));
                   objOut.d{k}.Z = objIn.d{k}.Z/sqrt(temp);
               end        
           else
               temp = 0;
               for k = 1:objIn.NTX
                   temp = temp+sum(sum(conj(objIn.d{k}.Z).*objIn.d{k}.Z));
               end  
               for k = 1:objIn.NTX
                   objOut.d{k}.Z = objIn.d{k}.Z/sqrt(temp);
               end        
           end
       end
       %**************************************************************************
       function obj = makeDataTemplate(obj,grid,xRows,yRows,T,type)
           %   make a data vector template for data of given type, with sites
           %   centered on cells defined by array rows x cols, for periods in T
           obj.NTX  = length(T);
           obj.d = cell(obj.NTX,1);
           %   data locations:
           
           x = [0  ; cumsum(grid.dx)]+grid.origin(1);
           y = [0 ; cumsum(grid.dy)]+grid.origin(2);
           xCtr = (x(1:end-1)+x(2:end))/2;
           yCtr = (y(1:end-1)+y(2:end))/2;
           xSite = xCtr(xRows);
           ySite = yCtr(yRows);
           [X,Y] = ndgrid(xSite,ySite);
           nSites = length(xSite)*length(ySite);
           siteLoc = [reshape(X,1,nSites);reshape(Y,1,nSites);zeros(1,nSites)];
           siteChar = num2str([1:nSites]');
           siteChar(siteChar==' ') = '0';
           
           switch upper(type)
               case 'FULL IMPEDANCE'
                   nComp = 8;
                   nTF = 4;
                   %compChar = ['Re(Zxx)';'Im(Zxx)';'Re(Zxy)';'Im(Zxy)'; ...
                   %   'Re(Zyx)';'Im(Zyx)';'Re(Zyy)';'Im(Zyy)'];
                   compChar = ['ZXX';'ZXY';'ZYX';'ZYY'];
               case 'OFF DIAGONAL'
                   nComp = 4;
                   nTF = 2;
                   compChar = ['ZXY';'ZYX'];
               case 'IMPEDANCE TIPPER'
                   nComp = 12;
                   nTF = 6;
                   compChar = ['ZXX';'ZXY';'ZYX';'ZYY';'TX ';'TY '];
                   %compChar = ['Re(Zxx)';'Im(Zxx)';'Re(Zxy)';'Im(Zxy)'; ...
                   %'Re(Zyx)';'Im(Zyx)';'Re(Zyy)';'Im(Zyy)'; ...
                   %'Re(Tx) ';'Im(Tx) ';'Re(Ty) ';'Im(Ty) '];
           end
           Z = zeros(nSites,nTF)+1i*zeros(nSites,nTF);
           Zerr = ones(nSites,nTF);
           
           for j=1:obj.NTX
               obj.d{j} = struct('T',T(j),...
                   'nComp',nComp,...
                   'compChar',compChar,...
                   'siteLoc',siteLoc',...
                   'siteChar',siteChar,...
                   'Z',Z,'Zerr',Zerr,...
                   'Cmplx',1, ...
                   'units','[V/m]/[T]',...
                   'signConventioni',-1,...
                    'origin',grid.origin,...
                    'orient',0,...
                    'lat',zeros(nSites,1),...
                    'lon',zeros(nSites,1));
           end
       end
   end    %    methods
end     %   classdef
