classdef MT3DmodelParam < ModelParameter
        % class for storing 2D MT model parameters

    properties %(SetAccess = private)
        v           %     array of conductive cells
        paramType   %   linear or log
        AirCond     %    air conductivity   
    end

    methods
       %*******************************************************************
       function obj = MT3DmodelParam(GRID,paramType)
       % class constructor
            if nargin == 1
                paramType='LOGE';
            end
            if nargin >=1
               nx = GRID.Nx;
               ny = GRID.Ny;
               nz = GRID.NzEarth;
               obj.v = zeros(nx,ny,nz);
               obj.paramType = paramType;
               obj.grid = GRID;
            end
       end
       %*******************************************************************
       function M = ModelParamLength(obj)
           %  sets MT2DmodelParam properties
           [nx,ny,nz] = size(obj.v);
           M = nx*ny*nz;
       end
       %*******************************************************************
       function M = ModelParamMax(obj)
           %  estimate maximum model parameter size ... used by NLCG for
           %  setting initial step size
           M = max(max(max(abs(obj.v))));
       end
       %*******************************************************************
       function obj = setModelParam(obj,v,paramType,AirCond,Grid)
           %  sets MT2DmodelParam properties
           obj.v = v;
           obj.paramType = paramType;
           obj.AirCond = AirCond;
           if nargin > 4
               obj.grid = Grid;
           end
       end
       %*******************************************************************
       function obj = setModAnomaly(obj,sigma0,anomalies)
           obj.v = sigma0*ones(size(obj.v));
           for k = 1:length(anomalies)
               i1 = anomalies{k}.xInd(1);
               i2 = anomalies{k}.xInd(2);
               j1 = anomalies{k}.yInd(1);
               j2 = anomalies{k}.yInd(2);
               k1 = anomalies{k}.zInd(1);
               k2 = anomalies{k}.zInd(2);
               obj.v(i1:i2,j1:j2,k1:k2) = anomalies{k}.sigma;
           end
           if strcmp(upper(obj.paramType),'LOGE')
               obj.v = log(obj.v);
           end
       end
       %*******************************************************************
       function grid = extractGrid(obj)
           grid = obj.grid;
       end
       %*******************************************************************
       function [obj] = readVec(obj,cfile)
       %   reads from ascii 2D model conductivity file cfile, puts
       %   contents into already created obj
       % Usage:  [m] = readCond(cfile)
       

       %    [dx,dy,dz,rho,nzAir,type,origin,rotation] = ...
       %        read_mackie3d_model(cfile);
           [dx,dy,dz,rho,nzAir,type,origin,rotation] = ...
                                             read_WS3d_model(cfile);
                                         
           if isempty(type)
               type = 'LINEAR';
           end

           Grid.dx = dx;
           Grid.dy = dy;
           Grid.dz = dz;
           Grid.Nx = length(dx);
           Grid.Ny = length(dy);
           Grid.NzEarth = length(dz);
           Grid.NzAir = nzAir;
           Grid.origin = origin;
           Grid.rotation = rotation;
           Grid.units = 'm';

           if strfind(type,'LOGE')
               obj.v = - rho;
           else
               obj.v = 1./rho;
           end
           obj.AirCond = log(1e-10);
           obj.paramType = type;
           obj.grid = Grid;

       end
       %*******************************************************************
       function [status] = writeVec(obj,cfile)
       %
       % Usage:  [status] = writeCond_2D(cfile,cond)
       %
       %  writes ModelParam object, provided as structure
       %  status is total number of bytes written
       
           %   always use Mackie format for conductivity parameters                 
           if ~isfield(obj.grid,'origin')
               obj.grid.origin = [0 0 0];
           end
           
           if isfield(obj.grid,'units')
               if strcmp(obj.grid.units,'km')
                   % convert everything to meters!
                   obj.grid.dx = 1000*obj.grid.dx;
                   obj.grid.dy = 1000*obj.grid.dy;
                   obj.grid.dz = 1000*obj.grid.dz;
                   obj.grid.origin = 1000*obj.grid.origin;
               end
           end
           if strcmp(obj.paramType,'LOGE')
               rho = - obj.v;
           else
               rho = 1./(obj.v);
           end
           
           status = write_WS3d_model(cfile,obj.grid.dx,obj.grid.dy,...
               obj.grid.dz,rho,obj.grid.NzAir,obj.paramType,...
               obj.grid.origin,obj.grid.rotation);
       end
       %*******************************************************************
       function mVec = ExtractVec(obj)
       %   extracts model parameter values and returns as a standard vector
           mVec = obj.v;
           [nx,ny,nz] = size(mVec);
           mVec = reshape(mVec,nx*ny*nz,1);
       end
        %*******************************************************************
       function obj = SetVec(obj,v)
       %   inserts  standard vector v into model parameter obj (as an
       %   array) 
           obj.v = ...
               reshape(v,[obj.grid.Nx,obj.grid.Ny,obj.grid.NzEarth]);
       end
       %*******************************************************************
       function obj = plus(obj1,obj2)
           obj = obj1;
           obj.v = obj1.v+obj2.v;
       end
       %*******************************************************************
       function obj = minus(obj1,obj2)
           obj = obj1;
           obj.v = obj1.v-obj2.v;
       end
       %*******************************************************************
       function [m] = linComb(c1,m1,c2,m2)
       %
       %   computes linear combination of model parameters
       %                 m = c1*m1+c2*m2
       %   This version is for a single model parameter
       %
       %  Usage: [m] = linCombMod(c1,m1,c2,m2);
           
           m = m1;
           m.v = c1*m1.v + c2*m2.v;
       end
       %*******************************************************************
       function [mOut] = mtimes(c,mIn)
           %  Initial implementation of scalar multiplication of a model
           %   space object by a real scalar; no error type checking;
           %   scalar has to be first argument
           %
           %  Usage: [mOut] = mtimes(c,mIn)
           if ~isa(c,'double')
               error('Error: can only multiply data vectors by doubles')
           end
           
           mOut = mIn;
           mOut.v = c*mOut.v;
       end
       %*******************************************************************
       function [ip] = dot(m1,m2)
       %  dot product for (real) model parameters m1, m2
       %
       %  Usage:  [ip] = dot(m1,m2);

          ip = sum(sum(sum(m1.v .* m2.v)));
       end
       %*******************************************************************
       function [mOut] = InitHalfSpace(mIn,sig)
       %   Usage: [mOut] = InitHalfSpace(mIn,sig);
       %   Initialize a model parameter mOut as half space
       %     mIn = template for model parameter
       %     sig = conductivity (NOT log conductivity,
       %		or resistivity) for half space
       %
       %  set prior, and starting conductivity (uniform half space)
           mOut = mIn;
           if strcmp(mIn.paramType,'LOGE')
               mOut.v = log(sig)*ones(size(mIn.v));
           else
               mOut.v = sig*ones(size(mIn.v));
           end
       end
       %*******************************************************************
       function [obj] = ZeroVec(obj)
       %   Usage: [mOut] = InitHalfSpace(mIn,sig);
       %    zero existing model parameter 
           obj.v = zeros(size(obj.v));
       end
       %*******************************************************************
       function plotCond(Cond)
       % Plots conductivity model
       % Usage [h,hCB] = plotCond(m)
       %   Cond is the 3D conductivity model parameter object to plot
       %    uses 3D_MT/COndPlot
      
           % convert to log10
           if strcmp(Cond.paramType,'LOGE')
               Cond.paramType = 'LOG10';
               Cond.v = Cond.v / log(10);
               Cond.AirCond = Cond.AirCond / log(10);
           elseif strcmp(Cond.paramType,'LINEAR')
               Cond.paramType = 'LOG10';
               Cond.v = log10(Cond.v);
               Cond.AirCond = log10(Cond.AirCond);
           end

           options.slice = 'Z';
           options.Np = 1;
           [Nx,Ny,Nz] = size(Cond.v);
           options.iXlim(1) = 1;
           options.iXlim(2) = Nx+1;
           options.iYlim(1) = 1;
           options.iYlim(2) = Ny+1;
           options.iZlim(1) = 1;
           options.iZlim(2) = Nz+1;
           options.cblabel = 'log_{10} \sigma';
           CondPlotSet(Cond.v,Cond.grid,options);
       
       end
       %***********************************************
       function diff  = relDiff(obj1,obj2)
           temp = obj1-obj2;
           diff = dot(temp,temp)/dot(obj1,obj1);
           diff = sqrt(diff);
       end
       %***********************************************
       function [d] = fwd(m,d0)
           %  Computes predicted data for conductivity sigma0
           %   Writes sigma0, dIn to files to be read by Fortran program Mod2DMT
           %   calls program with option to  compute predicted data, reads result
           %    and returns result in dOut
           %
           %  Usage: [d] = fwd(m,d0);
           %
           %  Inputs:	m = conductivity parameter (structure)
           %		d0 = data vector (cell array of structures)
           %  Output:	d = predicted data vector (cell array of structures)
           
           %    make scratch directory if it doesn't exist ...
           if exist('scratch','dir')==0
               mkdir('scratch');
           end
           
           %   write out conductivity parameter in scratch directory
           m_File = 'scratch/Input.cpr';
           writeVec(m,m_File);
           
           %   write out data vector template in scratch directory
           d0_File = 'scratch/Input.imp';
           writeVec(d0,d0_File);
           %   file name for computed data vector, to be output by Mod2DMT
           d_File = 'scratch/Out.imp';
           
           %   run Mod3DMT
           Test3D('FORWARD',m_File,d0_File,d_File)
           
           %  read in predicted data
           d = MT3DZ;
           d = readVec(d,d_File);
       end
    end    % methods
end   %   classdef
