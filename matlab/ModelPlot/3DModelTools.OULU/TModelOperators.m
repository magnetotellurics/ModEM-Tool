classdef TModelOperators < TGrid3D
% Egbert,Smirnov,Cherevatova 2013    
%
% defines differential operators of 3D rectangular grid
% Grad, Curl
% Fields are descritized on the edges

properties
  FaceArea             %   array defining face areas 
  NFaces
  EdgeLength             %   array defining edge lengths 
  NEdges
  NNodes
  DualFaceArea            %   array defining dual face areas (order as for E)
  DualEdgeLength            %   array defining dual edge lengths (order as for F)
  VNode           %   array of node volumes
  VCell           %   array of cell volumes 
  T             %   topology of curl operator: sparse matrix mapping
                %    edges to facesatheros
  G             %   topology of gradient operator: sparse matrix mapping
                %     nodes to edges
  DDx, DDy, DDz % dual grid cell dimensions
         
end
methods
function obj = TModelOperators(Dx,Dy,Dz,Nza)
  %  class constructor
  obj = obj@TGrid3D(Dx,Dy,Dz,Nza);
  obj.setCellVolume;  
  obj.setEdgeLength;
  obj.setFaceArea;
  obj.dualLengths;
  obj.setNodeVolume;
  obj.setDualEdgeLength;
  obj.setDualFaceArea;                    
  obj.setGrad;
  obj.setCurl;
end
%*******************************************************************        
function dualLengths(obj)
 %   compute dual grid length elements
  obj.DDx = ([obj.Dx; 0]+[0; obj.Dx])/2;
  obj.DDy = ([obj.Dy; 0]+[0; obj.Dy])/2;
  obj.DDz = ([obj.Dz; 0]+[0 ;obj.Dz])/2;
end

 %******************************************************************
function setCellVolume(obj)
  %   create 1-D array of cell volumes, ordered as indC
  [DX,DY,DZ] = ndgrid(obj.Dx,obj.Dy,obj.Dz);
  obj.VCell = DX.*DY.*DZ;
end
%
%******************************************************************
function setNodeVolume(obj)
  obj.NNodes = (obj.Nx+1)*(obj.Ny+1)*(obj.Nz+1);
  
 %   create 1-D array of node volumens, ordered as indN
  %    Let's just do this with a short-cut here, so that it
  %    is sensible for the special case where all cells in each
  %    vertical layer have the same coarseness
  %    SET edge lengths first

  [DX,DY,DZ] = ndgrid(obj.DDx,obj.DDy,obj.DDz);
  obj.VNode = DX.*DY.*DZ;            
 end
%
%******************************************************************
function setEdgeLength(obj)
  %   creates 1-D array of edge lengths
  nxe = obj.Nx*(obj.Ny+1)*(obj.Nz+1);
  nye = (obj.Nx+1)*obj.Ny*(obj.Nz+1);
  nze = (obj.Nx+1)*(obj.Ny+1)*obj.Nz;
  obj.NEdges = [nxe;nye;nze];
  obj.EdgeLength = zeros(nxe+nye+nze,1);
  %   x-edges   
  i1 = 1; i2 = obj.NEdges(1);
  obj.EdgeLength(i1:i2) = ...
      reshape(obj.Dx*ones(1,(obj.Ny+1)*(obj.Nz+1)),nxe,1);
  % y-edges
  i1 = i2+1; i2 = i2 + obj.NEdges(2);
  temp = reshape(obj.Dy*ones(1,obj.Nz+1),1,obj.Ny*(obj.Nz+1));
  obj.EdgeLength(i1:i2) = reshape(ones(obj.Nx+1,1)*temp,nye,1);
  %  z-edges
  i1 = i2+1; i2 = i2 + obj.NEdges(3);
  obj.EdgeLength(i1:i2) = reshape(ones((obj.Nx+1)*(obj.Ny+1),1)*obj.Dz',nze,1);  
end
%
%******************************************************************
function setFaceArea(obj)
  %   creates 1-D array of face areas
  nxf = (obj.Nx+1)*obj.Ny*(obj.Nz);
  nyf = obj.Nx*(obj.Ny+1)*(obj.Nz);
  nzf = obj.Nx*obj.Ny*(obj.Nz+1);

  obj.NFaces = [nxf;nyf;nzf];            
  obj.FaceArea = zeros(nxf+nyf+nzf,1);
  
  %   x-face elements
  i1 = 1; i2 = obj.NFaces(1);
  temp  = reshape(ones(obj.Nx+1,1)*obj.Dy',(obj.Nx+1)*obj.Ny,1);
  obj.FaceArea(i1:i2) = reshape(temp*obj.Dz',obj.Nx+1,obj.Ny,obj.Nz);

  %   y-face elements
  i1 = i2+1; i2 = i2+obj.NFaces(2);
  temp  = reshape(obj.Dx*ones(1,obj.Ny+1),obj.Nx*(obj.Ny+1),1);
  obj.FaceArea(i1:i2) = reshape(temp*obj.Dz',obj.Nx,obj.Ny+1,obj.Nz);          
  
  %  z-face elements
  i1 = i2+1; i2 = i2+obj.NFaces(3);
  temp = reshape(obj.Dx*obj.Dy',obj.Nx*obj.Ny,1);
  obj.FaceArea(i1:i2) = reshape(temp*ones(1,obj.Nz+1),obj.Nx,obj.Ny,obj.Nz+1);
end

 % sum(obj.nEd
 %******************************************************************
 function setDualEdgeLength(obj)
  %   creates 1-D array of dual edge lengths, one for each face         
  obj.DualEdgeLength = zeros(sum(obj.NFaces),1);
  %   x-edges 
  ne = obj.NFaces(1);
  i1 = 1; i2 = ne;
  obj.DualEdgeLength(i1:i2) = ...
      reshape(obj.DDx*ones(1,obj.Ny*obj.Nz),ne,1);
  % y-edges
  ne = obj.NFaces(2);
  i1 = i2+1; i2 = i2 + ne;
  temp = reshape(obj.DDy*ones(1,obj.Nz),1,(obj.Ny+1)*obj.Nz);
  obj.DualEdgeLength(i1:i2) = reshape(ones(obj.Nx,1)*temp,ne,1);
  %  z-edges
  ne = obj.NFaces(3);
  i1 = i2+1; i2 = i2 + ne;
  obj.DualEdgeLength(i1:i2) = reshape(ones(obj.Nx*obj.Ny,1)*obj.DDz',ne,1);  
end
 %
 %******************************************************************
function setDualFaceArea(obj)
   %   creates 1-D array of dual face areas, one for each edge
  obj.DualFaceArea = zeros(sum(obj.NEdges),1);
  
  %   x-face elements
  nf = obj.NEdges(1);
  i1 = 1; i2 = nf;
  temp  = reshape(ones(obj.Nx,1)*obj.DDy',obj.Nx*(obj.Ny+1),1);
  obj.DualFaceArea(i1:i2) = reshape(temp*obj.DDz',nf,1);

  %   y-face elements
  nf = obj.NEdges(2);
  i1 = i2+1; i2 = i2+nf;
  temp  = reshape(obj.DDx*ones(1,obj.Ny),(obj.Nx+1)*obj.Ny,1);
  obj.DualFaceArea(i1:i2) = reshape(temp*obj.DDz',nf,1);          
  
  %  z-face elements
  nf = obj.NEdges(3);
  i1 = i2+1; i2 = i2+nf;
  temp = reshape(obj.DDx*obj.DDy',(obj.Nx+1)*(obj.Ny+1),1);
  obj.DualFaceArea(i1:i2) = reshape(temp*ones(1,obj.Nz),nf,1);
             
end
%
%******************************************************************
function setGrad(obj)
    %   set sparse matrix defining gradient operator: maps from
    %   nodes to edges
    % x-edges
  n = obj.NEdges(1);
  i1 = 1; i2 = n;
  rowIndicies = reshape(ones(2,1)*(i1:i2),n*2,1);
  %   columns
  [I,J,K] = obj.gridIndex(1:n,'xedge');
  I2 = I+1;
  colIndicies = reshape( ...
    [obj.vectorIndex(I,J,K,'node')'; obj.vectorIndex(I2,J,K,'node')'], n*2,1);
  
  %   matrix entries
  entries = reshape([-ones(1,n); ones(1,n)],n*2,1);     
        
    % y-edges
  n = obj.NEdges(2);
  i1 = i2+1; i2 = i2+n;
  rowIndicies = [rowIndicies; reshape(ones(2,1)*(i1:i2),n*2,1)];
  %   columns
  [I,J,K] = obj.gridIndex(1:n,'yedge');
  J2 = J+1;
  colIndicies = [ colIndicies; ...
     reshape([ obj.vectorIndex(I,J,K,'node')'; obj.vectorIndex(I,J2,K,'node')'], n*2,1)];
  %   matrix entries
  entries = [entries; reshape([-ones(1,n); ones(1,n)],n*2,1)];
        
    % z-edges
  n = obj.NEdges(3);
  i1 = i2+1; i2 = i2+n;
  rowIndicies = [rowIndicies; reshape(ones(2,1)*(i1:i2),n*2,1)];
   %   columns
  [I,J,K] = obj.gridIndex(1:n,'zedge');
   %    quadtree: no subdivision of vertical layers
  K2 = K+1;
  colIndicies = [ colIndicies; 
    reshape( [obj.vectorIndex(I,J,K,'node')'; obj.vectorIndex(I,J,K2,'node')'], n*2,1)];
   %   matrix entries
  entries = [entries; reshape([-ones(1,n); ones(1,n)],n*2,1)];
        
%  define sparse matrix     
  obj.G = sparse(rowIndicies,colIndicies,entries);    
end %setGrad
%
%******************************************************************
function setCurl(obj)
   %   set sparse matrix defining curl operator: maps from
   %   edges to faces
  
   %   set sparse matrix defining curl operator: maps from
   %   edges to faces
   Nx = obj.Nx;
   Ny = obj.Ny;
   Nz = obj.Nz;
   
   i1 = 1; i2 = obj.NEdges(1);
   xEdgeTemp = zeros(Nx,Ny+1,Nz+1);
   xEdgeTemp(1:obj.NEdges(1)) = i1:i2;
   i1 = i2+1; i2 = i2+obj.NEdges(2);
   yEdgeTemp = zeros(Nx+1,Ny,Nz+1);
   yEdgeTemp(1:obj.NEdges(2)) = i1:i2;
   i1 = i2+1; i2 = i2+obj.NEdges(3);
   zEdgeTemp = zeros(Nx+1,Ny+1,Nz);
   zEdgeTemp(1:obj.NEdges(3)) = i1:i2;
            
%   x-faces
   n = obj.NFaces(1);
   i1 = 1; i2 = n;
   rowIndicies = reshape(ones(4,1)*(i1:i2),n*4,1);
   %   columns (for y-edge, z-edge)
   [I,J,K] = obj.gridIndex(1:n,'xface');
   J2 = J + 1;
   K2 = K + 1;
   colIndicies = reshape( ...
       [ yEdgeTemp(obj.vectorIndex(I,J,K,'yedge')); ...
         zEdgeTemp(obj.vectorIndex(I,J,K,'zedge')); ...
         yEdgeTemp(obj.vectorIndex(I,J,K2,'yedge')); ...
         zEdgeTemp(obj.vectorIndex(I,J2,K,'zedge'))], n*4,1);
   entries = [-ones(1,n) ; ones(1,n) ; ones(1,n) ; -ones(1,n)];
   entries = reshape(entries,4*n,1);
 

%   y-faces
   n = obj.NFaces(2);
   i1 = i2+1; i2 = i2+n;
   rowIndicies = [rowIndicies; ...
       reshape(ones(4,1)*(i1:i2),n*4,1)];
   %   columns (for y-edge, z-edge)
   [I,J,K] = obj.gridIndex(1:n,'yface');
   I2 = I+1;
   K2 = K+1;
   colIndicies = [colIndicies ; reshape( ...
       [ zEdgeTemp(obj.vectorIndex(I,J,K,'zedge')); ...
         xEdgeTemp(obj.vectorIndex(I,J,K,'xedge')); ...
         zEdgeTemp(obj.vectorIndex(I2,J,K,'zedge')); ...
         xEdgeTemp(obj.vectorIndex(I,J,K2,'xedge'))], n*4,1)];
   entries = [ entries ; ...
               reshape( [-ones(1,n) ; ones(1,n) ; ones(1,n) ; -ones(1,n)], n*4,1)];

   
%   z-faces
   n = obj.NFaces(3);
   i1 = i2+1; i2 = i2+n;
   rowIndicies = [rowIndicies; ...
       reshape(ones(4,1)*(i1:i2),n*4,1)];
   %   columns (for y-edge, z-edge)
   [I,J,K] = obj.gridIndex(1:n,'zface');
   I2 = I+1;
   J2 = J+1;
   colIndicies = [colIndicies ; reshape( ...
       [ xEdgeTemp(obj.vectorIndex(I,J,K,'xedge')); ...
         yEdgeTemp(obj.vectorIndex(I,J,K,'yedge')); ...       
         xEdgeTemp(obj.vectorIndex(I,J2,K,'xedge')); ...
         yEdgeTemp(obj.vectorIndex(I2,J,K,'yedge'))], n*4,1)];
   entries = [ entries ; ...
               reshape([-ones(1,n) ; ones(1,n) ; ones(1,n) ; -ones(1,n)], n*4,1)];
   
   ncol = length(i1:i2);
   nrow = length(i1:i2);
   obj.T = sparse(rowIndicies,colIndicies,entries);             
end % setCurl


% function setCurlCurl(obj)   
%     
%     
%  % all of the following is tested
%  % on non uniform grid
%  
% %  	Grad : E^{-1}*G
% %         Div  : V_d^{-1}*G'*F_d
% %         so Div * Grad ( = nabla!) V_d^{-1} * G' * F_d * E^{-1} * G
% %              This is symmetric after multiplication by the appropriate volume
% %                   elements (for nodes, as Div * Grad maps from nodes to nodes.
% % 	Curl = F^{-1} * T * E
% %         CurlT = F_d^{-1}*T'*E_d
% %         curl-curl  = CurlT * Curl = F_d^{-1} * T' * E_d * F^{-1} * T * E
% %             This is not what you have; as for Div * Grad this is symmetric only
% %               after being multiplied by the appropriate volume elements, this
% %                time the EDGE VOLUMES = V_e = E*F_d  (i.e., the volume of the cube
% %                 from cell center to cell center, for the four cells surrounding
% %                  the edge).
% % 
%  
%  
%  
%   % Volume weights  
%   V(1:obj.NNodes) = obj.VNode;  
%   
%   % Div
%   n = sum(obj.NEdges);
%  % obj.Div = spdiags(1./V', 0,obj.NNodes,obj.NNodes)*obj.G'*spdiags(obj.DualFaceArea,0,n,n);    
%   obj.VDiv = obj.G'*spdiags(obj.DualFaceArea,0,n,n);    
%   
%   % Grad  obj.E = sparse(obj.E);
% 
%   %obj.Grad = diag(1./obj.EdgeLength)*obj.G;
%   
%   
%   % non symmetric
%   %divCgrad = spdiags(1./V, 0,obj.NNodes,obj.NNodes)*obj.G'*diag(obj.DualFaceArea./obj.EdgeLength)*obj.G;
%   
%   %the following is symmetric by applying volume weitghs
%   %V_d*Div*Grad  =   G' * F_d * E^{-1} * G
%   %DivGrad = obj.G'*diag(obj.DualFaceArea./obj.EdgeLength)*obj.G;
%   
%   
%   %curl = diag(1./obj.FaceArea)*obj.T*diag(obj.EdgeLength);
%   
%   %curlT = diag(1./obj.DualFaceArea)*obj.T'*diag(obj.DualEdgeLength);
%   
%   % super symmetric
%   % symmetric after multiplying by volume weights
%   % V_e *curl*curl =  E * F_d *  F_d^{-1} * T' * E_d * F^{-1} * T * E 
%   %  =                E          * T'            E_d * F^{-1} *T *E 
%   ne = length(obj.EdgeLength);
%   nd = length(obj.DualEdgeLength);
%   obj.VCurlCurl = ...
%     spdiags(obj.EdgeLength,0,ne,ne)*obj.T'* ...
%     spdiags(obj.DualEdgeLength./obj.FaceArea,0,nd,nd)*obj.T*spdiags(obj.EdgeLength,0,ne,ne)  ;
%   
% end;
% 


end % methods
end % classdef
