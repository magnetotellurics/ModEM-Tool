classdef T3DfwdFDE < TModelOperators
% 
% object to direved from basic 
% model operators class
% extends to contain conductivities


properties
  omega     % W 2*pi/T
  sigmaCell  % array of cell conductivities
  sigmaEdge  % array of edge conductivities
  mu    % just single value for now  
  epsilon % could be an array too
  b     % RHS vector containing boundary conditions    
  A     % LHS
  ind_b % indeces of boundary edges
  ind_i % indeces of interior edges
  ind_n % indeces of inner nodes 
  Ex  % 3D array of Ez field
  Ey  % 3D array of Ez field
  Ez  % 3D array of Ez field
  E   % vector solution Ed
  Eb1D % vector of 1D boundary solution
  Eb2D % 2D boundary
  M   % sigma cell mapping to edges operator
  VDiv           % divergence operator V_d*Div
  Div
  VCurlCurl      % V_d*Curl'*Curl
  VDivCGrad
end

methods
%*******************************************************************
function obj = T3DfwdFDE(Dx,Dy,Dz,Nza)
  %   class constructor ... simple
  obj = obj@TModelOperators(Dx,Dy,Dz,Nza);  
  % allocate E components  
  obj.Ex = zeros(obj.Nx,obj.Ny+1,obj.Nz+1);  
  obj.Ey = zeros(obj.Nx+1,obj.Ny,obj.Nz+1); 
  obj.Ez = zeros(obj.Nx+1,obj.Ny+1,obj.Nz);  
  
% find ind of interior and boundary edges by setting
% boundaries to 1
  obj.Ex(1:obj.Nx,[1 obj.Ny+1],1:obj.Nz+1) = 1;  
  obj.Ex(1:obj.Nx,1:obj.Ny+1,[1 obj.Nz+1]) = 1;  
  
  obj.Ey([1 obj.Nx+1],1:obj.Ny,1:obj.Nz+1) = 1;
  obj.Ey(1:obj.Nx+1,1:obj.Ny,[1 obj.Nz+1]) = 1;
  
  obj.Ez(1:obj.Nx+1,[1 obj.Ny+1],1:obj.Nz) = 1;  
  obj.Ez([1 obj.Nx+1],1:obj.Ny+1,1:obj.Nz) = 1;  

  % create solution vector E 
  obj.mapE2vec;
  % boundary edges
  obj.ind_b = find(obj.E);
  % inner edges
  obj.ind_i = find(obj.E == 0);
  % inner nodes indices
  i = reshape(repmat((2:obj.Nx)', obj.Ny-1,obj.Nz-1)  ,(obj.Ny-1)*(obj.Nz-1)*(obj.Nx-1),1);
  j = reshape(repmat((2:obj.Ny), obj.Ny-1,obj.Nz-1)   ,(obj.Ny-1)*(obj.Nz-1)*(obj.Nx-1),1);
  k = reshape(ones((obj.Ny-1)*(obj.Nx-1),1)*(2:obj.Nz),(obj.Ny-1)*(obj.Nz-1)*(obj.Nx-1),1);
  
  obj.ind_n = obj.vectorIndex(i,j,k,'node')';
  % set default
  obj.mu = pi*4E-7;
  obj.epsilon = 1;
  obj.omega = 2*pi*1000;   % 1000Hz
end
%*******************************************************************

function setSigmaMap(obj)
  % used for mapping of the model parameters from edges to cells   
  % set mapping matrix M from edge to cell
  % unweighted cell2edge is
  % SigmaCell =  SigmaEdge(ix, iy, iz)/4.0
  % setMapping = 1 1 1 1 ...divided by 4
   Nx = obj.Nx;
   Ny = obj.Ny;
   Nz = obj.Nz;
   NCells = Nx*Ny*Nz;
   
   i1 = 1; i2 = obj.NEdges(1);
   xEdgeTemp = zeros(Nx,Ny+1,Nz+1);
   xEdgeTemp(1:obj.NEdges(1)) = i1:i2;
   i1 = i2+1; i2 = i2+obj.NEdges(2);
   yEdgeTemp = zeros(Nx+1,Ny,Nz+1);
   yEdgeTemp(1:obj.NEdges(2)) = i1:i2;
   i1 = i2+1; i2 = i2+obj.NEdges(3);
   zEdgeTemp = zeros(Nx+1,Ny+1,Nz);
   zEdgeTemp(1:obj.NEdges(3)) = i1:i2;
 
  % x-edges
  n = NCells;
  i1 = 1; i2 = n;
  rowIndicies = reshape(ones(4,1)*(i1:i2),n*4,1);
  %   columns
  [I,J,K] = obj.gridIndex(1:n,'cell');
  J2 = J + 1;
  K2 = K + 1;

  colIndicies = reshape( ...
     [xEdgeTemp(obj.vectorIndex(I,J,K,'xedge')); ...
      xEdgeTemp(obj.vectorIndex(I,J,K2,'xedge')); ...
      xEdgeTemp(obj.vectorIndex(I,J2,K,'xedge')); ...
      xEdgeTemp(obj.vectorIndex(I,J2,K2,'xedge'))], n*4,1);


   entries = [ones(1,n)/4. ; ones(1,n)/4. ; ones(1,n)/4. ; ones(1,n)/4.];
   entries = reshape(entries,4*n,1);

  % y-edges
  %i1 = i2+1; i2 = i2+n;
  rowIndicies =  [rowIndicies; reshape(ones(4,1)*(i1:i2),n*4,1)];
  %   columns
  [I,J,K] = obj.gridIndex(1:n,'cell');
  I2 = I + 1;
  K2 = K + 1;

  colIndicies = [colIndicies ;reshape( ...
    [yEdgeTemp(obj.vectorIndex(I,J,K,'yedge')); ...
     yEdgeTemp(obj.vectorIndex(I,J,K2,'yedge')); ...
     yEdgeTemp(obj.vectorIndex(I2,J,K,'yedge')); ...
     yEdgeTemp(obj.vectorIndex(I2,J,K2,'yedge'))], n*4,1)];

entries = [ entries ; ...
            reshape( [ones(1,n)/4. ; ones(1,n)/4. ; ones(1,n)/4. ; ones(1,n)/4.], n*4,1)];

  % z-edges
 % i1 = i2+1; i2 = i2+n;
  rowIndicies =  [rowIndicies; reshape(ones(4,1)*(i1:i2),n*4,1)];
  %   columns
  [I,J,K] = obj.gridIndex(1:n,'cell');
  I2 = I + 1;
  J2 = J + 1;
  colIndicies = [colIndicies ;reshape( ...
    [zEdgeTemp(obj.vectorIndex(I,J,K,'zedge')); ...
     zEdgeTemp(obj.vectorIndex(I,J2,K,'zedge')); ...
     zEdgeTemp(obj.vectorIndex(I2,J,K,'zedge')); ...
     zEdgeTemp(obj.vectorIndex(I2,J2,K,'zedge'))], n*4,1)];
 
    entries = [ entries ; ...
            reshape( [ones(1,n)/4. ; ones(1,n)/4. ; ones(1,n)/4. ; ones(1,n)/4.], n*4,1)];
      
   ncol = length(i1:i2);
   nrow = length(i1:i2);
   obj.M = sparse(colIndicies,rowIndicies,entries);             
     
end% setSigmaMap(obj)

function mapSigma2Edge(obj)
  % maps model parameters (Sigma or log(Sigma))
  % from cell centers to edges, where fields are defined
  % for symmetry the edge conductivity should respresent
  % volume weighted average
  % Sigma_edge = V_edge^(-1) * M * V_cell *Sigma_cell
  % where Sigma_edge -- edge conductivity
  %       V_edge^(-1) -- edge volumes
  %       M -- unweighted mapping matrix from edges to cells
  %                   setSigmaMap(obj)
  %       V_cell -- cell volumes
  %       Sigma_cell -- cell conductivity 
  % acts on both inner and boundary nodes
 
  ne =  length(obj.EdgeLength);
  nv =  obj.Nx*obj.Ny*obj.Nz;
  TempsigmaCell = reshape(obj.sigmaCell, obj.Nx*obj.Ny*obj.Nz,1);
  TempVcell = reshape(obj.VCell,obj.Nx*obj.Ny*obj.Nz,1);   
  
  % volume average
  obj.sigmaEdge = obj.M*(TempVcell.*TempsigmaCell)./(obj.M*TempVcell);

  %simple average 
  %obj.sigmaEdge = obj.M*TempsigmaCell;
        
end;    

function LoadModelFromFile(obj, fname)  
% reads a 3D resistivity model in Weerachai Siripunvaraporn's "0" format
if nargin == 2

 fid = fopen(fname);
 line = fgetl(fid); % comment
 line = fgetl(fid);
 [n] = sscanf(line,'%d',[4 1]);
 if findstr(line,'LOGE')
     type = 'LOGE';
 else
     type = 'LINEAR';
 end
 obj.Nx =   n(1);   obj.Ny  =   n(2);   obj.Nz  =   n(3);
 obj.Nza = 10;
 obj.Dx   =   fscanf(fid,'%f',obj.Nx);
 obj.Dy   =   fscanf(fid,'%f',obj.Ny);
 obj.Dz   =   fscanf(fid,'%f',obj.Nz);
 obj.Dz   =   [ zeros(obj.Nza,1); obj.Dz];
 k = 1;
 rho(1:obj.Nx,1:obj.Ny,1:obj.Nz+obj.Nza) = 0;
 while k <= obj.Nz
     for j = 1:obj.Ny
         for i = obj.Nx:-1:1
             obj.sigmaCell(i,j,k+obj.Nza) = 1./fscanf(fid,'%f',1);
         end
     end
     k = k+1;
 end
 line = fgetl(fid); % read to the end of line
 fclose(fid);
 %set Air 
 obj.sigmaCell(:,:,1:obj.Nza) = 1E-8;
 obj.Nz = obj.Nz + obj.Nza;
 obj.Dz(1:obj.Nza) = obj.Dz(obj.Nza+1)*2.^(obj.Nza-1:-1:0);
% if obj.Dz(1) < 30000; obj.Dz(1) = 30000; end;

% origin = [-sum(x)/2 -sum(y)/2 0];
% rotation = 0;
% while 1
%      line = fgetl(fid);
%      if ~ischar(line); 
%          fclose(fid);
%          return; 
%      end  
%     [n] = sscanf(line,'%f',[3 1]);
%     if length(n)==3
%         origin = n;
%     else
%         [n] = sscanf(line,'%f',1);
%         if length(n)==1
%             rotation = n; 
%         end
%     end    
%     
% here we can load Sigma from file        
  
else
 obj.sigmaEdge(1:sum(obj.NEdges),1) = 1;
 obj.sigmaCell(1:obj.Nx, 1:obj.Ny, 1:obj.Nz) = 1;
end; 
  
  obj.mu = 4*pi*1E-7;
  obj.epsilon = 1;
  obj.omega = 2*pi*1000;    
end;  

function setEquations(obj)   
    
    
% all of the following is tested
% on non uniform grid 
%  	Grad : E^{-1}*G
%         Div  : V_d^{-1}*G'*F_d
%         so Div * Grad ( = nabla!) V_d^{-1} * G' * F_d * E^{-1} * G
%              This is symmetric after multiplication by the appropriate volume
%                   elements (for nodes, as Div * Grad maps from nodes to nodes.
% 	Curl = F^{-1} * T * E
%         CurlT = F_d^{-1}*T'*E_d
%         curl-curl  = CurlT * Curl = F_d^{-1} * T' * E_d * F^{-1} * T * E
%             This is not what you have; as for Div * Grad this is symmetric only
%               after being multiplied by the appropriate volume elements, this
%                time the EDGE VOLUMES = V_e = E*F_d  (i.e., the volume of the cube
%                 from cell center to cell center, for the four cells surrounding
%                  the edge).
%  
 
  % Div
%  n = sum(obj.NEdges);
 % obj.VDiv = obj.G'*spdiags(obj.DualFaceArea,0,n,n);    
  
  %Grad  obj.E = sparse(obj.E);
  %obj.Grad = diag(1./obj.EdgeLength)*obj.G;
   
  % non symmetric
  %divCgrad = spdiags(1./V, 0,obj.NNodes,obj.NNodes)*obj.G'*diag(obj.DualFaceArea./obj.EdgeLength)*obj.G;
  
  %the following is symmetric by applying volume weitghs
  %V_d*Div*Grad  =   G' * F_d * E^{-1} * G
  %DivGrad = obj.G'*diag(obj.DualFaceArea./obj.EdgeLength)*obj.G;
    
  %curl = diag(1./obj.FaceArea)*obj.T*diag(obj.EdgeLength);
  
  %curlT = diag(1./obj.DualFaceArea)*obj.T'*diag(obj.DualEdgeLength);
  
  % super symmetric
  % symmetric after multiplying by volume weights
  % V_e *curl*curl =  E * F_d *  F_d^{-1} * T' * E_d * F^{-1} * T * E 
  %  =                E          * T'            E_d * F^{-1} *T *E 
  ne = length(obj.EdgeLength);
  nd = length(obj.DualEdgeLength);
  obj.VCurlCurl = ...
    spdiags(obj.EdgeLength,0,ne,ne)*obj.T'* ...
    spdiags(obj.DualEdgeLength./obj.FaceArea,0,nd,nd)*obj.T*spdiags(obj.EdgeLength,0,ne,ne)  ;
  
  % symmetric system V_d (Curl Curl inner + ioms) = V_d b       
  Vioms = -1i*obj.omega*obj.mu*obj.sigmaEdge(obj.ind_i).* ...
          obj.EdgeLength(obj.ind_i).*obj.DualFaceArea(obj.ind_i); 
  n = length(Vioms);
  obj.A = obj.VCurlCurl(obj.ind_i, obj.ind_i) + spdiags(Vioms,0,n,n);    


  % equations for divergence correction
  % 
  
  n = length(obj.DualFaceArea(obj.ind_i));
   
  obj.VDiv = obj.G(obj.ind_i,obj.ind_n)'*spdiags(obj.DualFaceArea(obj.ind_i).* ...
  obj.sigmaEdge(obj.ind_i), 0 ,n,n);

  obj.VDivCGrad = obj.G(obj.ind_i,obj.ind_n)'*spdiags(obj.DualFaceArea(obj.ind_i)./obj.EdgeLength(obj.ind_i).*obj.sigmaEdge(obj.ind_i), 0 ,n,n)* ...
  obj.G(obj.ind_i,obj.ind_n);  

  nn = length(obj.VNode(obj.ind_n));
  obj.Div = spdiags(1./reshape(obj.VNode(obj.ind_n), nn,1), 0,nn,nn)* ...
      obj.G(obj.ind_i,obj.ind_n)'*spdiags(obj.DualFaceArea(obj.ind_i).*obj.sigmaEdge(obj.ind_i), 0 ,n,n);
    
end;

function  mapE2vec(obj)
  % creates 1-D array of E fields
  % the ordering is exactly the same as for Edges
  obj.E = [reshape(obj.Ex,obj.NEdges(1),1); ...
           reshape(obj.Ey,obj.NEdges(2),1); ...
           reshape(obj.Ez,obj.NEdges(3),1)];
end;


function  mapE2grid(obj)
  % Ex
  i1 = 1; i2 = obj.NEdges(1);
  obj.Ex = reshape(full(obj.E(i1:i2)), obj.Nx,(obj.Ny+1),(obj.Nz+1));
  % Ey
  i1 = i2+1; i2 = i2 + obj.NEdges(2);
  obj.Ey = reshape(full(obj.E(i1:i2)), obj.Nx+1, obj.Ny, obj.Nz+1);
  %  Ez
  i1 = i2+1; i2 = i2 + obj.NEdges(3);
  obj.Ez = reshape(full(obj.E(i1:i2)), obj.Nx+1, obj.Ny+1, obj.Nz);
end;


function set1DBoundary(obj, Polarization)
% b is already mulpiplied by V_d volume weights
% V_d*b
  %z  = cumsum(obj.Dz(obj.Nza+1:obj.Nz)); 
  %k = (1+1i)*sqrt((obj.mu*obj.sigmaCell(1,1,obj.Nza+1)*obj.omega)/2);
  %obj.Eb1D = [ones(obj.Nza+1,1); exp(-k*z)];
  Nz = obj.Nz;
  sigmaCell1D = squeeze(obj.sigmaCell(1,1,1:Nz));
  Dz = obj.Dz;
  G = ( diag(-ones(Nz,1)) + diag( ones(Nz-1,1),1) );
  A = G'*diag(Dz.^-2)*G - 1i*obj.omega*obj.mu*diag(sigmaCell1D);

  b = -A(2:Nz,1);
  A = sparse(A(2:Nz,2:Nz));
  [L,U] = ilu(A,struct('type','ilutp','droptol',1e-6));  
  [x,flag,relres2,iter2,resvec2] = qmr(A,b,1e-20,100, L,U);
  obj.Eb1D = [1; x];
  obj.Eb1D(Nz+1) = 1E-10;
 
  obj.Ez = zeros(obj.Nx+1,obj.Ny+1,obj.Nz);  
  obj.Ex = zeros(obj.Nx,obj.Ny+1,obj.Nz+1);  
  obj.Ey = zeros(obj.Nx+1,obj.Ny,obj.Nz+1); 
  switch Polarization
      case 'X'
        for i=1:obj.Ny+1  
          obj.Ex(:,i,:) = repmat(obj.Eb1D',obj.Nx,1);
        end;
      case 'Y'              
        for i=1:obj.Nx+1  
          obj.Ey(i,:,:) = repmat(obj.Eb1D',obj.Ny,1);
        end;      
  end;
  obj.mapE2vec; 
  obj.b = -obj.VCurlCurl(obj.ind_i,obj.ind_b)*obj.E(obj.ind_b);
  
  
  
      
end;



function result = DivCorrection(obj)
    
  
  b = obj.VDiv*obj.E(obj.ind_i);
  A = obj.VDivCGrad;
  

  [x,flag,relres2,iter2,resvec2] = pcg(A,b,1e-6,100);
  deltaE =  diag(1./obj.EdgeLength(obj.ind_i))*obj.G(obj.ind_i,obj.ind_n)*x;
  obj.E(obj.ind_i) = obj.E(obj.ind_i) - deltaE;
  
  
  result = norm(obj.Div*obj.E(obj.ind_i));  
end;


function Solve3D(obj)
    
 maxiter = 1; iter = 0; notconverged = 1; tolerance = 1E-6;  
% [L,U] = ilu(obj.A, struct('type','ilutp','droptol',1e-6) );  
 [L,U] = ilu(obj.A);
 while (iter<maxiter) && notconverged   
   Eit = obj.E; 
   [x,flag,relres2,iter2,resvec2] = qmr(obj.A,obj.b,1e-20,10, L,U, obj.E(obj.ind_i));
%   [x,flag,relres2,iter2,resvec2] = qmr(obj.A,obj.b,1e-20,10, L,U);
%   [x,flag,relres2,iter2,resvec2] = qmr(obj.A,obj.b,1e-6,10);
   obj.E(obj.ind_i) = x;
   obj.mapE2grid;
 %  d = obj.DivCorrection;
 %  obj.mapE2grid; 
   d =0;
   nn = norm(Eit-obj.E);
   notconverged = norm(Eit-obj.E) > tolerance;
   iter=iter+1;
   fprintf('iter  %d  norm %f  div %f \n', iter,nn,d);
 end;
 
end;


end     % methods
end    % classdef