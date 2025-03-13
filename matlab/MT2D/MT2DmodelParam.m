classdef MT2DmodelParam < ModelParameter
        % class for storing 2D MT model parameters; a specific instance of
        % the abstract ModelParameter superclass

    properties %(SetAccess = private)
        v           %     array of conductive cells
        paramType   %   linear or log
        AirCond     %    air conductivity   
    end

    methods
       %*******************************************************************
       function obj = MT2DmodelParam(GRID,paramType)
       % class constructor
            if nargin == 1
                paramType='LOGE';
            end
            if nargin >=1
               ny = GRID.Ny;
               nz = GRID.Nz-GRID.Nza;
               obj.v = zeros(ny,nz);
               obj.paramType = paramType;
               obj.grid = GRID;
            end
       end
       %*******************************************************************
       function M = ModelParamLength(obj)
           %  sets MT2DmodelParam properties
           [ny,nz] = size(obj.v);
           M = ny*nz;
       end
       %*******************************************************************
       function M = ModelParamMax(obj)
           %  estimate maximum model parameter size ... used by NLCG for
           %  setting initial step size
           M = max(max(abs(obj.v)));
       end
       %*******************************************************************
       function obj = setModelParam(obj,v,paramType,AirCond)
           %  sets MT2DmodelParam properties
           obj.v = v;
           obj.paramType = paramType;
           obj.AirCond = AirCond;
       end
       %*******************************************************************
       function obj = SetFrom3D(obj,Cond,XY,slice)     
       %  fill in the contents of a 2D MT model parameter object from a
       %       3D conductivity model  input structure
           gr = MT2Dgrid;
           MT2DgridFrom3Dgrid(gr,Cond.grid,XY)
           obj.grid = gr;
           obj.paramType = 'LOGE';
           if strcmp(Cond.paramType,'LOG10')
               obj.AirCond = log(10)*Cond.AirCond;
               if strcmpi(XY,'x')
                   obj.v = log(10)*squeeze(Cond.v(slice,:,:));
               else
                   obj.v = log(10)*squeeze(Cond.v(:,slice,:));   
               end
           else
               obj.AirCond = log(Cond.AirCond);
               if strcmpi(XY,'x')
                   obj.v = log(squeeze(Cond.v(slice,:,:)));
               else
                   obj.v = log(squeeze(Cond.v(:,slice,:)));   
               end
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
       

           [dy,dz,rho,type] = read_mackie2d_model(cfile);

           if isempty(type)
               type = 'LINEAR';
           end

           ny = length(dy);
           nzEarth = length(dz);
           nzAir = 10;

           dzAir(1:nzAir) = 0;
           dzAir(1) = max(dz(1),10.0);
           for k = 2:nzAir
               dzAir(k) = dzAir(k-1)*3;
           end

           Grid.Dy = dy';
           Grid.Dz = [dzAir(end:-1:1) dz'];
           Grid.Ny = ny;
           Grid.Nz = nzEarth + nzAir;
           Grid.Nza = nzAir;

           if findstr(type,'LOGE')
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
           
           nza = obj.grid.Nza;
           nz  = obj.grid.Nz;
           
           dy = obj.grid.Dy;
           dz = obj.grid.Dz(nza+1:nz);
           
           type = obj.paramType;
           
           if strcmp(type,'LOGE')
               rho = - obj.v;
           else
               rho = 1./(obj.v);
           end
           
           status = write_mackie2d_model(cfile,dy,dz,rho,type);
       end
       %*******************************************************************
       function mVec = ExtractVec(obj)
       %   extracts model parameter values and returns as a standard vector
           mVec = obj.v;
           [ny,nz] = size(mVec);
           mVec = reshape(mVec,ny*nz,1);
       end
        %*******************************************************************
       function obj = SetVec(obj,v)
       %   inserts  standard vector v into model parameter obj (as an
       %   array) 
           v = reshape(v,[obj.grid.Ny,obj.grid.Nz-obj.grid.Nza]);
           obj.v = v;
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

          ip = sum(sum(m1.v .* m2.v));
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
       %   Usage: [mOut] = ZeroVec(obj);
       %    zero existing model parameter 
           obj.v = zeros(size(obj.v));
       end
       %*******************************************************************
       function [mOut] = oneD_average(mIn)
       %   Usage: [mOut] = oneD_averae(mIn);
       %    average mIn to a 1-D model
           mOut = mIn;
           mOut.v = ones(mIn.grid.Ny,1)*mean(mIn.v,1);
       end
       %*******************************************************************
       function [d] = fwd(m,d0)
           %  Computes predicted data for conductivity sigma0
           %   Writes sigma0, dIn to files to be read by Fortran program Mod2DMT
           %   calls program with option to  compute predicted data, reads result
           %    and returns result in dOut
           %
           %  Usage: [d] = MT2Dfwd(m,d0);
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
           
           %   run Mod2DMT
           Test2D('FORWARD',m_File,d0_File,d_File)
           
           %  read in predicted data
           d = MT2DZ;
           d = readVec(d,d_File);
       end
       %*******************************************************************
       function [h,hCB] = plotCond(m,OPTIONS)
       % Plots conductivity model
       % Usage [h,hCB] = plotCond(m,grid,OPTIONS)
       %   m is the 2D conductivity model to plot
       %   OPTIONS is a structure of plotting OPTIONS
       %          .nySkip = number of cells to omit from
       %                    each end of the profile
       %          .nZplot = number of vertical layers to plot
       %                    (starting from the earth surface)
       %          .ncax = color axis limits (vector with 2 elements)
       %          .title = plot title
           
           if nargin ==  1
               nZplot = m.grid.Nz-m.grid.Nza;
               nYskip = 5;
               cax = [-3.,0];
               OPTIONS = struct('Contour',0,'nYskip',nYskip,'nZplot',nZplot,...
                   'cax',cax,'title','');
           end
           if ~isfield(OPTIONS,'Contour')
               OPTIONS.Contour = 0;
           end
           if strcmp(m.paramType,'LOGE')
               m.v = log10(exp(m.v));
           else
               m.v = log10(m.v);
           end
           
           figure('Position',[100,100,600,400],...
               'PaperPosition',[1,1,6,4])
           y = cumsum([0 m.grid.Dy])';
           z = cumsum([0 m.grid.Dz])';
           zCenter = z-z(m.grid.Nza+1);
           zCenter = zCenter(m.grid.Nza+1:end)/1000;
           yCenter = y(OPTIONS.nYskip+1:end-OPTIONS.nYskip)/1000;
           yCenter = yCenter-mean(yCenter);
           Cond = m.v;
           Cond = [Cond Cond(:,end)];
           Cond = [Cond ;Cond(end,:)];
           Cond = Cond(OPTIONS.nYskip+1:end-OPTIONS.nYskip,1:OPTIONS.nZplot);
           if OPTIONS.Contour
               %   positive contours: solid line
               cL  = OPTIONS.cLev;
               cL = cL(cL>0); 
               [h,hCB] = contour(yCenter,zCenter(1:OPTIONS.nZplot),Cond',...
                   cL,'Linewidth',1,'color','b');
               set(gca,'FontWeight','bold','FontSize',16);
               clabel(h,hCB,cL(2:2:end),'FontSize',14,'fontweight','bold');
               hold on
               cL  = OPTIONS.cLev;
               cL = cL(cL<0);
               [h,hCB] = contour(yCenter,zCenter(1:OPTIONS.nZplot),Cond',...
                   cL,'Linewidth',1,'color','r','linestyle','--');
               clabel(h,hCB,cL(2:2:end),'FontSize',14,'fontweight','bold');
               [h,hCB] = contour(yCenter,zCenter(1:OPTIONS.nZplot),Cond',...
                   [0,0],'Linewidth',2,'Color',[0,0,0]);
               if length(h)>0
                  clabel(h,hCB,'manual','FontSize',14,'fontweight','bold')
               end
               axis('ij');
           else
              h = pcolor(yCenter,zCenter(1:OPTIONS.nZplot),Cond'); ...
                  shading flat; 
              axis('ij');
              set(gca,'FontWeight','demi','FontSize',13);
              caxis(OPTIONS.cax);
              c = colmap;
              c = c(end:-1:1,:);
              X = 1:17;
              XI = 1:.25:17;
              c = interp1(X,c,XI);
              colormap(c)
              hCB = colorbar;
              yt = floor(OPTIONS.cax(1)):1:ceil(OPTIONS.cax(2));
              ytLabel = 10.^(-yt);
              
              set(hCB,'FontWeight','demi','FontSize',12,...
                  'Ytick',yt,'YtickLabel',ytLabel);
           end
           ylabel('Depth (km)');
           xlabel('km');
           title(OPTIONS.title);
       end
       %***********************************************
       function diff  = relDiff(obj1,obj2)
           temp = obj1-obj2;
           diff = dot(temp,temp)/dot(obj1,obj1);
           diff = sqrt(diff);
       end
    end    % methods
end   %   classdef
