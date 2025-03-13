classdef quadtreeGrid < handle
    properties
        coarseGrid    %  handle of grid3D object; defines basic coarse grid
        %    individual vertical layers can be subdivided
        %    into quadtree grids
        Kmax             %  maximum depth of quadtree grid: max number of
        %    subdivisions (of a side) is 2^Kmax
        fineGrid      %   handle of grid3D object; defines the maximally
        %    subdivided fine grid
        Cells         %  3D array (fine grid dimensions) defining size
        indC          %   index of non-zero cells in full array of cells
        %    (in powers of 2) of defined cells, undefined
        %    cells as zeros
        CellsSet = false
        Nodes         %  As for Cells; (Nx+1,Ny+1,Nz+1)
        indN          %   index of non-zero nodes in full array of nodes
        NodesSet = false
        xFaces        %  (Nx+1,Ny,Nz)
        yFaces        %  (Nx,Ny+1,Nz)      These are all big arrays with lots
        zFaces        %  (Nx,Ny,Nz+1)         of zeros (potentially)
        indF          % index of non-zero elements in full array
                      %    of faces    [xFaces; yFaces ; zFaces] 
        nFaces
        FacesSet = false
        xEdges        %  (Nx,Ny+1,Nz+1)     Not all of these need to be
        yEdges        %  (Nx+1,Ny,Nz+1)        defined/carried around
        zEdges        %  (Nx+1,Ny+1,Nz)
        indE          %  index of non-zero elements in full array
                      %    of edges    [xEdges; yEdges ; zEdges]
        nEdges
        EdgesSet = false
        %   not sure if the actual operators should be defined here
        F             %   array defining face areas   
        E             %   array defining edge lengths  
        Fd            %   array defining dual face areas (order as for E)
        Ed            %   array defining dual edge lengths (order as for F)
        V_N           %   array of node volumes
        V_C           %   array of cell volumes 
        T             %   topology of curl operator: sparse matrix mapping
                      %    edges to faces, using orderings in indE and indF
        D             %   topology of gradient operator: sparse matrix mapping
                      %     nodes to edges, using orderings in indN, indE
    end
    methods
        function obj = quadtreeGrid(coarseGrid,Kmax)
            %  class constructor
            if nargin >= 1
                obj.coarseGrid = coarseGrid;
                obj.Kmax = 0;
            end
            if nargin == 2
                obj.Kmax = Kmax;
                obj.fineGrid = refineGridQT(obj.coarseGrid,Kmax);
            end
        end
        %
        %******************************************************************
        function setCellsLayered(obj,Klayer)
            %   given quadtreeGrid object, and grid depth for each vertical
            %   layer, construct Cells array
            if abs(diff(Klayer))>1
                error('adjacent layers can only change depth by one')
            end
            if max(Klayer)~=obj.Kmax
                k = max(Klayer);
                obj.Kmax = k;
                obj.fineGrid = refineGridQT(obj.coarseGrid,k);
            end
            obj.Cells = zeros(obj.fineGrid.Nx,obj.fineGrid.Ny,...
                obj.fineGrid.Nz);
            for k = 0:obj.Kmax
                indZ = Klayer == k;
                kk = obj.Kmax-k;
                indX = 1:2^kk:obj.fineGrid.Nx;
                indY = 1:2^kk:obj.fineGrid.Ny;
                obj.Cells(indX,indY,indZ) = 2^kk;
                %   fill in all fine cells aside from upper left with
                %   -2^k ... allows us to tell from any subcell the degree
                %    of coarseness of the actual cell containing this 
                %    fine-grid cell.
                for i=1:2^kk
                    for j=1:2^kk
                        if i~=1 || j~=1
                        obj.Cells(indX+i-1,indY+j-1,indZ) = -2^kk;
                        end
                    end
                end
            end
            obj.CellsSet = true;
            obj.indC = find(obj.Cells>0);
        end
        %
        %******************************************************************
        function setNodes(obj)
            %  set Node array, using Cells array
            if ~obj.CellsSet
                error('Cant set Nodes until Cells are set')
            end
            obj.Nodes = zeros(obj.fineGrid.Nx+1,obj.fineGrid.Ny+1,...
                obj.fineGrid.Nz+1);
            Nz = obj.fineGrid.Nz;
            for k = obj.Kmax:-1:0
                indX = 1:2^k:obj.fineGrid.Nx;
                indY = 1:2^k:obj.fineGrid.Ny;
                Nxk = length(indX);
                Nyk = length(indY);
                Nodes_k = zeros(Nxk+1,Nyk+1,Nz+1);
                % indicator array of cells that are at depth k
                kInd = obj.Cells(indX,indY,:) == 2^k;
                kInd = reshape(kInd,Nxk,Nyk,obj.fineGrid.Nz);
                Nodes_k(1:end-1,1:end-1,1:end-1) = kInd;
                Nodes_k(2:end,1:end-1,1:end-1) = ...
                    Nodes_k(2:end,1:end-1,1:end-1)+kInd;
                Nodes_k(1:end-1,2:end,1:end-1) = ...
                    Nodes_k(1:end-1,2:end,1:end-1)+kInd;
                Nodes_k(1:end-1,1:end-1,2:end) = ...
                    Nodes_k(1:end-1,1:end-1,2:end)+kInd;
                Nodes_k(2:end,2:end,1:end-1) = ...
                    Nodes_k(2:end,2:end,1:end-1)+kInd;
                Nodes_k(2:end,1:end-1,2:end) = ...
                    Nodes_k(2:end,1:end-1,2:end)+kInd;
                Nodes_k(1:end-1,2:end,2:end) = ...
                    Nodes_k(1:end-1,2:end,2:end)+kInd;
                Nodes_k(2:end,2:end,2:end) = ...
                    Nodes_k(2:end,2:end,2:end)+kInd;
                indX = 1:2^k:obj.fineGrid.Nx+1;
                indY = 1:2^k:obj.fineGrid.Ny+1;
                temp = obj.Nodes(indX,indY,:);
                temp(Nodes_k>0) = 0;
                obj.Nodes(indX,indY,:) = temp+2^k*(Nodes_k > 0);
            end
            obj.NodesSet = true;
            obj.indN = find(obj.Nodes>0);
        end
        %
        %******************************************************************
        function setEdges(obj)
            %  set Edge arrays, using Cells array
            if ~obj.CellsSet
                error('Cant set Edges until Cells are set')
            end
            obj.xEdges = zeros(obj.fineGrid.Nx,obj.fineGrid.Ny+1,...
                obj.fineGrid.Nz+1);
            obj.yEdges = zeros(obj.fineGrid.Nx+1,obj.fineGrid.Ny,...
                obj.fineGrid.Nz+1);
            obj.zEdges = zeros(obj.fineGrid.Nx+1,obj.fineGrid.Ny+1,...
                obj.fineGrid.Nz);
            Nz = obj.fineGrid.Nz;
            for k = obj.Kmax:-1:0
                indX = 1:2^k:obj.fineGrid.Nx;
                indY = 1:2^k:obj.fineGrid.Ny;
                Nxk = length(indX);
                Nyk = length(indY);
                xEdges_k = zeros(Nxk,Nyk+1,Nz+1);
                yEdges_k = zeros(Nxk+1,Nyk,Nz+1);
                zEdges_k = zeros(Nxk+1,Nyk+1,Nz);
                kInd = obj.Cells(indX,indY,:) == 2^k;
                kInd = reshape(kInd,Nxk,Nyk,obj.fineGrid.Nz);
                
                xEdges_k(1:end,1:end-1,1:end-1) = kInd;
                xEdges_k(1:end,2:end,1:end-1) = ...
                    xEdges_k(1:end,2:end,1:end-1)+kInd;
                xEdges_k(1:end,1:end-1,2:end) = ...
                    xEdges_k(1:end,1:end-1,2:end)+kInd;
                xEdges_k(1:end,2:end,2:end) = ...
                    xEdges_k(1:end,2:end,2:end)+kInd;
                indX = 1:2^k:obj.fineGrid.Nx;
                indY = 1:2^k:obj.fineGrid.Ny+1;
                temp = obj.xEdges(indX,indY,:);
                temp(xEdges_k>0) = 0;
                obj.xEdges(indX,indY,:) = temp+2^k*(xEdges_k > 0);
                
                yEdges_k(1:end-1,1:end,1:end-1) = kInd;
                yEdges_k(2:end,1:end,1:end-1) = ...
                    yEdges_k(2:end,1:end,1:end-1)+kInd;
                yEdges_k(1:end-1,1:end,2:end) = ...
                    yEdges_k(1:end-1,1:end,2:end)+kInd;
                yEdges_k(2:end,1:end,2:end) = ...
                    yEdges_k(2:end,1:end,2:end)+kInd;
                indX = 1:2^k:obj.fineGrid.Nx+1;
                indY = 1:2^k:obj.fineGrid.Ny;
                temp = obj.yEdges(indX,indY,:);
                temp(yEdges_k>0) = 0;
                obj.yEdges(indX,indY,:) = temp+2^k*(yEdges_k > 0);
                
                zEdges_k(1:end-1,1:end-1,1:end) = kInd;
                zEdges_k(2:end,1:end-1,1:end) = ...
                    zEdges_k(2:end,1:end-1,1:end)+kInd;
                zEdges_k(1:end-1,2:end,1:end) = ...
                    zEdges_k(1:end-1,2:end,1:end)+kInd;
                zEdges_k(2:end,2:end,1:end) = ...
                    zEdges_k(2:end,2:end,1:end)+kInd;
                indX = 1:2^k:obj.fineGrid.Nx+1;
                indY = 1:2^k:obj.fineGrid.Ny+1;
                temp = obj.zEdges(indX,indY,:);
                temp(zEdges_k>0) = 0;
                %   unclear how z-edges (which are really all the full
                %     coarse layer thickness in this quadtree scheme (i.e.,
                %     no subdivision of layers)   should these all be 1?
                obj.zEdges(indX,indY,:) = temp+2^k*(zEdges_k > 0);        
            end
            obj.nEdges = zeros(3,1);
            obj.EdgesSet = true;
            temp = find(obj.xEdges>0);
            obj.nEdges(1) = length(temp);
            obj.indE = temp;
            temp = find(obj.yEdges>0);
            obj.nEdges(2) = length(temp);
            obj.indE = [obj.indE ; temp];
            temp = find(obj.zEdges>0);
            obj.nEdges(3) = length(temp);
            obj.indE = [obj.indE ; temp];
        end
        %
        %******************************************************************
        function setFaces(obj)
            %  set Face arrays, using Cells array
            if ~obj.CellsSet
                error('Cant set Faces until Cells are set')
            end
            obj.xFaces = zeros(obj.fineGrid.Nx+1,obj.fineGrid.Ny,...
                obj.fineGrid.Nz);
            obj.yFaces = zeros(obj.fineGrid.Nx,obj.fineGrid.Ny+1,...
                obj.fineGrid.Nz);
            obj.zFaces = zeros(obj.fineGrid.Nx,obj.fineGrid.Ny,...
                obj.fineGrid.Nz+1);
            %   x-faces
            Nz = obj.fineGrid.Nz;
            for k = obj.Kmax:-1:0
                indX = 1:2^k:obj.fineGrid.Nx;
                indY = 1:2^k:obj.fineGrid.Ny;
                Nxk = length(indX);
                Nyk = length(indY);
                xFaces_k = zeros(Nxk+1,Nyk,Nz);
                yFaces_k = zeros(Nxk,Nyk+1,Nz);
                zFaces_k = zeros(Nxk,Nyk,Nz+1);
                kInd = obj.Cells(indX,indY,:) == 2^k;
                kInd = reshape(kInd,Nxk,Nyk,obj.fineGrid.Nz);
                
                xFaces_k(1:end-1,1:end,1:end) = kInd;
                xFaces_k(2:end,1:end,1:end) = ...
                    xFaces_k(2:end,1:end,1:end)+kInd;
                indX = 1:2^k:obj.fineGrid.Nx+1;
                indY = 1:2^k:obj.fineGrid.Ny;
                temp = obj.xFaces(indX,indY,:);
                temp(xFaces_k>0) = 0;
                obj.xFaces(indX,indY,:) = temp+2^k*(xFaces_k > 0);
                
                yFaces_k(1:end,1:end-1,1:end) = kInd;
                yFaces_k(1:end,2:end,1:end) = ...
                    yFaces_k(1:end,2:end,1:end)+kInd;
                indX = 1:2^k:obj.fineGrid.Nx;
                indY = 1:2^k:obj.fineGrid.Ny+1;
                temp = obj.yFaces(indX,indY,:);
                temp(yFaces_k>0) = 0;
                obj.yFaces(indX,indY,:) = temp+2^k*(yFaces_k > 0);
                
                zFaces_k(1:end,1:end,1:end-1) = kInd;
                zFaces_k(1:end,1:end,2:end) = ...
                    zFaces_k(1:end,1:end,2:end)+kInd;
                indX = 1:2^k:obj.fineGrid.Nx;
                indY = 1:2^k:obj.fineGrid.Ny;
                temp = obj.zFaces(indX,indY,:);
                temp(zFaces_k>0) = 0;
                obj.zFaces(indX,indY,:) = temp+2^k*(zFaces_k > 0);
            end
            obj.FacesSet = true;
            obj.nFaces = zeros(3,1);
            temp = find(obj.xFaces>0);
            obj.nFaces(1) = length(temp);
            obj.indF = temp;
            temp = find(obj.yFaces>0);
            obj.nFaces(2) = length(temp);
            obj.indF = [obj.indF ; temp];
            temp = find(obj.zFaces>0);
            obj.nFaces(3) = length(temp);
            obj.indF = [obj.indF ; temp];
        end
        %
        %******************************************************************
        function setCellVolume(obj)
            %   create 1-D array of cell volumens, ordered as indC
            [DX,DY,DZ] = ndgrid(obj.fineGrid.Dx,obj.fineGrid.Dy,obj.fineGrid.Dz);
            temp = (abs(obj.Cells).^2).*DX.*DY.*DZ;
            obj.V_C = temp(obj.indC);
        end
        %
        %******************************************************************
        function setEdgeLength(obj)
            %   creates 1-D array of edge lengths, orderd as indE
            Nx = obj.fineGrid.Nx;
            Ny = obj.fineGrid.Ny;
            Nz = obj.fineGrid.Nz;
            %   x-edges
            temp = reshape(obj.fineGrid.Dx*ones(1,(Ny+1)*(Nz+1)),Nx,Ny+1,Nz+1);
            temp = temp.*obj.xEdges;
            obj.E = zeros(length(obj.indE),1);
            i1 = 1; i2 = obj.nEdges(1);
            obj.E(i1:i2) = temp(obj.indE(i1:i2));
            %  y-edges
            temp = reshape(obj.fineGrid.Dy*ones(1,Nz+1),1,Ny*(Nz+1));
            temp = reshape(ones(Nx+1,1)*temp,Nx+1,Ny,Nz+1);
            temp = temp.*obj.yEdges;
            i1 = i2+1; i2 = i2+obj.nEdges(2);
            obj.E(i1:i2) = temp(obj.indE(i1:i2));
            %  z-edges
            temp = reshape(ones((Nx+1)*(Ny+1),1)*obj.fineGrid.Dz',Nx+1,Ny+1,Nz);
            temp = temp.*(obj.zEdges>0);
            i1 = i2+1; i2 = i2+obj.nEdges(3);
            obj.E(i1:i2) = temp(obj.indE(i1:i2));
        end
        %
        %******************************************************************
        function setFaceArea(obj)
             %   creates 1-D array of face areas, orderd as indF
            Nx = obj.fineGrid.Nx;
            Ny = obj.fineGrid.Ny;
            Nz = obj.fineGrid.Nz;
            Dx = obj.fineGrid.Dx;
            Dy = obj.fineGrid.Dy;
            Dz = obj.fineGrid.Dz;
            obj.F = zeros(length(obj.indF),1);
            %   x-face elements
            temp  = reshape(ones(Nx+1,1)*Dy',(Nx+1)*Ny,1);
            temp = reshape(temp*Dz',Nx+1,Ny,Nz);
            %    Quadtree: z-edges are not subdivied
            temp = temp.*obj.xFaces;
            i1 = 1; i2 = obj.nFaces(1);
            obj.F(i1:i2) = temp(obj.indF(i1:i2));
            %   y-face elements
            temp  = reshape(Dx*ones(1,Ny+1),Nx*(Ny+1),1);
            temp = reshape(temp*Dz',Nx,Ny+1,Nz);
             %    Quadtree: z-edges are not subdivied
            temp = temp.*obj.yFaces;
            i1 = i2+1; i2 = i2+obj.nFaces(2);
            obj.F(i1:i2) = temp(obj.indF(i1:i2));
            %  z-face elements
            temp = reshape(Dx*Dy',Nx*Ny,1);
            temp = reshape(temp*ones(1,Nz+1),Nx,Ny,Nz+1);
            temp = temp.*obj.zFaces.^2;
            i1 = i2+1; i2 = i2+obj.nFaces(3);
            obj.F(i1:i2) = temp(obj.indF(i1:i2));
        end
        %
        %******************************************************************
        function setNodeVolume(obj)
            %   create 1-D array of node volumens, ordered as indN
            %     WRONG   WRONG WRONG   WRONG
            [delX,delY,delZ] = dualLengths(obj.fineGrid);
            [DX,DY,DZ] = ndgrid(delX,delY,delZ);
            temp = (obj.Nodes.^2)*DX.*DY.*DZ;
            obj.V_N = temp(obj.indN);
        end
        %
        %******************************************************************
        function setDualEdgeLength(obj)
            %   creates 1-D array of dual edge lengths, orderd as indF
            %    definition of dual grid is not clear for general case!!!
            
            Nx = obj.fineGrid.Nx;
            Ny = obj.fineGrid.Ny;
            Nz = obj.fineGrid.Nz;
            [dx,dy,dz] = ndgrid(obj.fineGrid.Dx,obj.fineGrid.Dy,obj.fineGrid.Dz);
            obj.Ed = zeros(length(obj.indF),1);
            
            %   dual x-edges
            temp = zeros(Nx+1,Ny,Nz);
            temp(1:end-1,:,:) = dx.*abs(obj.Cells);
            temp(2:end,:,:) = temp(2:end,:,:)+ dx.*abs(obj.Cells);
            temp = temp.*(obj.xFaces>0)/2;
            i1 = 1; i2 = obj.nFaces(1);
            obj.Ed(i1:i2) = temp(obj.indF(i1:i2));
            
            %   dual y-edges
            temp = zeros(Nx,Ny+1,Nz);
            temp(:,1:end-1,:) = dy.*abs(obj.Cells);
            temp(:,2:end,:) = temp(:,2:end,:)+ dy.*abs(obj.Cells);
            temp = temp.*(obj.yFaces>0)/2;
            i1 = 1+i2; i2 = i2+obj.nFaces(2);
            obj.Ed(i1:i2) = temp(obj.indF(i1:i2));
            
            %   dual z-edges
            temp = zeros(Nx,Ny,Nz+1);
            temp(:,:,1:end-1) = dz.*abs(obj.Cells);
            temp(:,:,2:end) = temp(:,:,2:end)+ dz.*abs(obj.Cells);
            temp = temp.*(obj.zFaces>0)/2;
            i1 = 1+i2; i2 = i2+obj.nFaces(3);
            obj.Ed(i1:i2) = temp(obj.indF(i1:i2));  
        end
        %
        %******************************************************************
        function setDualFaceArea(obj)
             %   creates 1-D array of face areas, orderd as indE
             %   Call after setDualEdgeLength -- uses these to compute area
             %    definition of dual grid is not clear for general case!!!
             
            Nx = obj.fineGrid.Nx;
            Ny = obj.fineGrid.Ny;
            Nz = obj.fineGrid.Nz;
            
            %    make temporary 3D arrays containing indicies into face
            %    vector; pad two edges with zeros to simplify calculation
            %    on edges
            i1 = 1; i2 = obj.nFaces(1);
            xFaceTemp = zeros(Nx+1,Ny+2,Nz+2);
            temp = zeros(Nx+1,Ny,Nz);
            temp(obj.indF(i1:i2)) = i1:i2;
            xFaceTemp(:,2:end-1,2:end-1) = temp;
            
            i1 = i2+1; i2 = i2+obj.nFaces(2);
            yFaceTemp = zeros(Nx+2,Ny+1,Nz+2);
            temp = zeros(Nx,Ny+1,Nz);
            temp(obj.indF(i1:i2)) = i1:i2;
            yFaceTemp(2:end-1,:,2:end-1) = temp;
            
            i1 = i2+1; i2 = i2+obj.nFaces(3);
            zFaceTemp = zeros(Nx+2,Ny+2,Nz+1);
            temp = zeros(Nx,Ny,Nz+1);
            temp(obj.indF(i1:i2)) = i1:i2;
            zFaceTemp(2:end-1,2:end-1,:) = temp;
            
            obj.Fd = zeros(length(obj.indE),1);
            
            %   dual x-faces
            %     find all faces bounded by each x-edge
            n = obj.nEdges(1);
            i1 = 1; i2 = n;
            [I,J,K] = gridIndex(obj.fineGrid,obj.indE(i1:i2),'xEdge');
            %   get dual edges associated with adjacent y-faces
                indDown = yFaceTemp(vectorIndex(obj.fineGrid,I+1,J+1,K+1,'yface'));
                dyDown = zeros(length(indDown),1);
                dyDown(indDown>0) = obj.Ed(indDown(indDown>0));
                %  next lines are to allow for the case where the face below is
                %  coarser than the edge, and the edge is the second
                ind = find(indDown==0);
                indDown2 = yFaceTemp(vectorIndex( ...
                    obj.fineGrid,I(ind),J(ind)+1,K(ind)+1,'yface'));
                dyDown(ind(indDown2>0)) = obj.Ed(indDown2(indDown2>0));

                indUp = yFaceTemp(vectorIndex(obj.fineGrid,I+1,J+1,K,'yface'));
                dyUp = zeros(length(indUp),1);
                dyUp(indUp>0) = obj.Ed(indUp(indUp>0));
                %  next lines are to allow for the case where the face above is
                %  coarser than the edge, and the edge is the second
                ind = find(indUp==0);
                indUp2 = yFaceTemp(vectorIndex( ...
                    obj.fineGrid,I(ind),J(ind)+1,K(ind),'yface'));
                dyUp(ind(indUp2>0)) = obj.Ed(indUp2(indUp2>0));
                dy = (dyUp+dyDown)/2;
                
            %   get dual edges associated with adjacent z-faces
                indLeft = zFaceTemp(vectorIndex(obj.fineGrid,I+1,J+1,K+1,'zface'));
                dzLeft = zeros(length(indLeft),1);
                dzLeft(indLeft>0) = obj.Ed(indLeft(indLeft>0));
                %  next lines are to allow for the case where the face to Left is
                %  coarser than the edge, and the edge is the second
                ind = find(indLeft==0);
                indLeft2 = zFaceTemp(vectorIndex( ...
                    obj.fineGrid,I(ind),J(ind)+1,K(ind)+1,'zface'));
                dzLeft(ind(indLeft2>0)) = obj.Ed(indLeft2(indLeft2>0));

                indRight = zFaceTemp(vectorIndex(obj.fineGrid,I+1,J,K+1,'zface'));
                dzRight = zeros(length(indRight),1);
                dzRight(indRight>0) = obj.Ed(indRight(indRight>0));
                %  next lines are to allow for the case where the face to Right is
                %  coarser than the edge, and the edge is the second
                ind = find(indRight==0);
                indRight2 = zFaceTemp(vectorIndex( ...
                    obj.fineGrid,I(ind),J(ind),K(ind)+1,'zface'));
                dzRight(ind(indRight2>0)) = obj.Ed(indRight2(indRight2>0));
                dz = (dzLeft+dzRight)/2;
                obj.Fd(i1:i2) = dy.*dz;
                
            %   dual y-faces
            %     find all faces bounded by each y-edge
            n = obj.nEdges(2);
            i1 = i2+1; i2 = i2+n;
            [I,J,K] = gridIndex(obj.fineGrid,obj.indE(i1:i2),'yEdge');
            %   get dual edges associated with adjacent x-faces
                indDown = xFaceTemp(vectorIndex(obj.fineGrid,I+1,J+1,K+1,'xface'));
                dxDown = zeros(length(indDown),1);
                dxDown(indDown>0) = obj.Ed(indDown(indDown>0));
                %  next lines are to allow for the case where the face below is
                %  coarser than the edge, and the edge is the second
                ind = find(indDown==0);
                indDown2 = xFaceTemp(vectorIndex( ...
                    obj.fineGrid,I(ind)+1,J(ind),K(ind)+1,'xface'));
                dxDown(ind(indDown2>0)) = obj.Ed(indDown2(indDown2>0));

                indUp = xFaceTemp(vectorIndex(obj.fineGrid,I+1,J+1,K,'xface'));
                dxUp = zeros(length(indUp),1);
                dxUp(indUp>0) = obj.Ed(indUp(indUp>0));
                %  next lines are to allow for the case where the face above is
                %  coarser than the edge, and the edge is the second
                ind = find(indUp==0);
                indUp2 = xFaceTemp(vectorIndex( ...
                    obj.fineGrid,I(ind)+1,J(ind),K(ind),'yface'));
                dxUp(ind(indUp2>0)) = obj.Ed(indUp2(indUp2>0));
                dx = (dxUp+dxDown)/2;
                
            %   get dual edges associated with adjacent z-faces
                indLeft = zFaceTemp(vectorIndex(obj.fineGrid,I,J+1,K+1,'zface'));
                dzLeft = zeros(length(indLeft),1);
                dzLeft(indLeft>0) = obj.Ed(indLeft(indLeft>0));
                %  next lines are to allow for the case where the face to Left is
                %  coarser than the edge, and the edge is the second
                ind = find(indLeft==0);
                indLeft2 = zFaceTemp(vectorIndex( ...
                    obj.fineGrid,I(ind),J(ind),K(ind)+1,'zface'));
                dzLeft(ind(indLeft2>0)) = obj.Ed(indLeft2(indLeft2>0));

                indRight = zFaceTemp(vectorIndex(obj.fineGrid,I+1,J+1,K+1,'zface'));
                dzRight = zeros(length(indRight),1);
                dzRight(indRight>0) = obj.Ed(indRight(indRight>0));
                %  next lines are to allow for the case where the face to Right is
                %  coarser than the edge, and the edge is the second
                ind = find(indRight==0);
                indRight2 = zFaceTemp(vectorIndex( ...
                    obj.fineGrid,I(ind)+1,J(ind),K(ind)+1,'zface'));
                dzRight(ind(indRight2>0)) = obj.Ed(indRight2(indRight2>0));
                dz = (dzLeft+dzRight)/2;
                obj.Fd(i1:i2) = dx.*dz;
 
            %   dual z-faces
            %     find all faces bounded by each z-edge ... 
            %        simpler case, as long as z levels are not subdivided
            n = obj.nEdges(3);
            i1 = i2+1; i2 = i2+n;
            [I,J,K] = gridIndex(obj.fineGrid,obj.indE(i1:i2),'zEdge');
            %   get dual edges associated with adjacent x-faces
                indRight = xFaceTemp(vectorIndex(obj.fineGrid,I+1,J,K+1,'xface'));
                dxRight = zeros(length(indRight),1);
                dxRight(indRight>0) = obj.Ed(indRight(indRight>0));
                indLeft = xFaceTemp(vectorIndex(obj.fineGrid,I+1,J+1,K+1,'xface'));
                dxLeft = zeros(length(indLeft),1);
                dxLeft(indLeft>0) = obj.Ed(indLeft(indLeft>0));
                dx = (dxRight+dxLeft)/2;
            %   get dual edges associated with adjacent y-faces 
                indRight = yFaceTemp(vectorIndex(obj.fineGrid,I+1,J+1,K+1,'yface'));
                dyRight = zeros(length(indRight),1);
                dyRight(indRight>0) = obj.Ed(indRight(indRight>0));
                indLeft = xFaceTemp(vectorIndex(obj.fineGrid,I,J+1,K+1,'yface'));
                dyLeft = zeros(length(indLeft),1);
                dyLeft(indLeft>0) = obj.Ed(indLeft(indLeft>0));
                dy = (dyRight+dyLeft)/2;
                obj.Fd(i1:i2) = dx.*dy;
        end
        %
        %******************************************************************
        function setDiv(obj)
            %   set sparse matrix defining divergence operator: maps from
            %   nodes to edges
            
            Nx = obj.fineGrid.Nx;
            Ny = obj.fineGrid.Ny;
            Nz = obj.fineGrid.Nz;
            %  need node indices stored in full array indexed by i,j,k
            NodeTemp = zeros(Nx+1,Ny+1,Nz+1);
            NodeTemp(obj.indN) = 1:length(obj.indN);
            % x-edges
                n = obj.nEdges(1);
                i1 = 1; i2 = n;
                rowIndicies = reshape(ones(2,1)*(i1:i2),n*2,1);
                %   columns
                [I,J,K] = gridIndex(obj.fineGrid,obj.indE(i1:i2),'xedge');
                I2 = I+obj.xEdges(obj.indE(i1:i2));
                colIndicies = reshape( ...
                    [NodeTemp(vectorIndex(obj.fineGrid,I,J,K,'node'))'; ...
                     NodeTemp(vectorIndex(obj.fineGrid,I2,J,K,'node'))'], ...
                      n*2,1);
                
                %   matrix entries
                entries = reshape([-ones(1,n); ones(1,n)],n*2,1);
            
            % y-edges
                n = obj.nEdges(2);
                i1 = i2+1; i2 = i2+n;
                rowIndicies = [rowIndicies; reshape(ones(2,1)*(i1:i2),n*2,1)];
                %   columns
                [I,J,K] = gridIndex(obj.fineGrid,obj.indE(i1:i2),'yedge');
                J2 = J+obj.yEdges(obj.indE(i1:i2));
                colIndicies = [ colIndicies; reshape( ...
                      [ NodeTemp(vectorIndex(obj.fineGrid,I,J,K,'node'))'; ...
                      NodeTemp(vectorIndex(obj.fineGrid,I,J2,K,'node'))'], ...
                      n*2,1)];
                %   matrix entries
                entries = [entries; reshape([-ones(1,n); ones(1,n)],n*2,1)];
            
            % z-edges
                n = obj.nEdges(3);
                i1 = i2+1; i2 = i2+n;
                rowIndicies = [rowIndicies; reshape(ones(2,1)*(i1:i2),n*2,1)];
                %   columns
                [I,J,K] = gridIndex(obj.fineGrid,obj.indE(i1:i2),'zedge');
                %    quadtree: no subdivision of vertical layers
                K2 = K+1;
                colIndicies = [ colIndicies; reshape( ...
                      [NodeTemp(vectorIndex(obj.fineGrid,I,J,K,'node'))'; ...
                      NodeTemp(vectorIndex(obj.fineGrid,I,J,K2,'node'))'], ...
                      n*2,1)];
                %   matrix entries
                entries = [entries; reshape([-ones(1,n); ones(1,n)],n*2,1)];
                
            
            %  define sparse matrix
             nrow = length(obj.indE);
             ncol = length(obj.indN);
            obj.D = sparse(rowIndicies,colIndicies,entries,nrow,ncol);  
        end
        %
        %******************************************************************
        function setCurl(obj)
            %   set sparse matrix defining curl operator: maps from
            %   edges to faces
            Nx = obj.fineGrid.Nx;
            Ny = obj.fineGrid.Ny;
            Nz = obj.fineGrid.Nz;
            i1 = 1; i2 = obj.nEdges(1);
            xEdgeTemp = zeros(Nx,Ny+1,Nz+1);
            xEdgeTemp(obj.indE(i1:i2)) = i1:i2;
            i1 = i2+1; i2 = i2+obj.nEdges(2);
            yEdgeTemp = zeros(Nx+1,Ny,Nz+1);
            yEdgeTemp(obj.indE(i1:i2)) = i1:i2;
            i1 = i2+1; i2 = i2+obj.nEdges(3);
            zEdgeTemp = zeros(Nx+1,Ny+1,Nz);
            zEdgeTemp(obj.indE(i1:i2)) = i1:i2;
            %   x-faces
                n = obj.nFaces(1);
                i1 = 1; i2 = n;
                rowIndicies = reshape(ones(4,1)*(i1:i2),n*4,1);
                %   columns (for y-edge, z-edge)
                [I,J,K] = gridIndex(obj.fineGrid,obj.indF(i1:i2),'xface');
                J2 = J+obj.xFaces(obj.indF(i1:i2));
                %    quadtree: no subdivision of vertical layers
                K2 = K+1;
                colIndicies = reshape( ...
                    [ yEdgeTemp(vectorIndex(obj.fineGrid,I,J,K,'yedge'))'; ...
                      zEdgeTemp(vectorIndex(obj.fineGrid,I,J,K,'zedge'))'; ...
                      yEdgeTemp(vectorIndex(obj.fineGrid,I,J,K2,'yedge'))'; ...
                      zEdgeTemp(vectorIndex(obj.fineGrid,I,J2,K,'zedge'))'], ...
                      n*4,1);
                entries = [-ones(1,n) ; ones(1,n) ; ones(1,n) ; -ones(1,n)];
                entries = reshape(entries,4*n,1);
                %   some additional y-edges need to be included, to allow
                %   for subdivided edges bounding a larger cell
                ind = find(obj.xFaces(obj.indF(i1:i2))> ...
                    obj.yEdges(vectorIndex(obj.fineGrid,I,J,K,'yedge')));
                n1 = length(ind);
                rows = ind;
                shift = obj.xFaces(obj.indF(i1:i2));
                I2  = I(ind)+shift(ind)/2;
                cols = yEdgeTemp(vectorIndex(obj.fineGrid,I2,J(ind),K(ind),'yedge'));
                entries = [entries ; -ones(n1,1)];
                rowIndicies = [rowIndicies; rows];
                colIndicies = [colIndicies; cols];
                %    y-edge along bottom of x-face
                ind = find(obj.xFaces(obj.indF(i1:i2))> ...
                    obj.yEdges(vectorIndex(obj.fineGrid,I,J,K2,'yedge')));
                n1 = length(ind);
                rows = ind;
                I2  = I(ind)+shift(ind)/2;
                cols = yEdgeTemp(vectorIndex(obj.fineGrid,I2,J(ind),K2(ind),'yedge'));
                entries = [entries ; ones(n1,1)];
                rowIndicies = [rowIndicies; rows];
                colIndicies = [colIndicies; cols];
                
             %   y-faces
                n = obj.nFaces(2);
                i1 = i2+1; i2 = i2+n;
                rowIndicies = [rowIndicies; ...
                    reshape(ones(4,1)*(i1:i2),n*4,1)];
                %   columns (for y-edge, z-edge)
                [I,J,K] = gridIndex(obj.fineGrid,obj.indF(i1:i2),'yface');
                I2 = I+obj.yFaces(obj.indF(i1:i2));
                %    quadtree: no subdivision of vertical layers
                K2 = K+1;
                colIndicies = [colIndicies ; reshape( ...
                    [ zEdgeTemp(vectorIndex(obj.fineGrid,I,J,K,'zedge'))'; ...
                      xEdgeTemp(vectorIndex(obj.fineGrid,I,J,K,'xedge'))'; ...
                      zEdgeTemp(vectorIndex(obj.fineGrid,I2,J,K,'zedge'))'; ...
                      xEdgeTemp(vectorIndex(obj.fineGrid,I,J,K2,'xedge'))'], ...
                      n*4,1)];
                entries = [ entries ; reshape( ...
                    [-ones(1,n) ; ones(1,n) ; ones(1,n) ; -ones(1,n)], ...
                    n*4,1)];
                %   some additional x-edges need to be included, to allow
                %   for subdivided edges bounding a larger cell
                
                ind = find(obj.yFaces(obj.indF(i1:i2)) > ...
                     obj.xEdges(vectorIndex(obj.fineGrid,I,J,K,'xedge')));
                n1 = length(ind);
                rows = ind + obj.nFaces(1);
                shift = obj.yFaces(obj.indF(i1:i2));
                J2  = J(ind)+shift(ind)/2;
                cols = xEdgeTemp(vectorIndex(obj.fineGrid,I(ind),J2,K(ind),'xedge'));
                entries = [entries ; ones(n1,1)];
                rowIndicies = [rowIndicies; rows];
                colIndicies = [colIndicies; cols];
                %    x-edge along bottom of x-face
                ind = find(obj.yFaces(obj.indF(i1:i2)) > ...
                     obj.xEdges(vectorIndex(obj.fineGrid,I,J,K2,'xedge')));
                n1 = length(ind);
                rows =ind + obj.nFaces(1);
                J2  = J(ind)+shift(ind)/2;
                cols = xEdgeTemp(vectorIndex(obj.fineGrid,I(ind),J2,K2(ind),'xedge'));
                entries = [entries ; -ones(n1,1)];
                rowIndicies = [rowIndicies; rows];
                colIndicies = [colIndicies; cols];
                
                %   z-faces
                n = obj.nFaces(3);
                i1 = i2+1; i2 = i2+n;
                rowIndicies = [rowIndicies; ...
                    reshape(ones(4,1)*(i1:i2),n*4,1)];
                %   columns (for y-edge, z-edge)
                [I,J,K] = gridIndex(obj.fineGrid,obj.indF(i1:i2),'zface');
                I2 = I+obj.zFaces(obj.indF(i1:i2));
                J2 = J+obj.zFaces(obj.indF(i1:i2));
                colIndicies = [colIndicies ; reshape( ...
                    [ xEdgeTemp(vectorIndex(obj.fineGrid,I,J,K,'xedge'))'; ...
                    yEdgeTemp(vectorIndex(obj.fineGrid,I,J,K,'yedge'))'; ...
                    
                    xEdgeTemp(vectorIndex(obj.fineGrid,I,J2,K,'xedge'))'; ...
                    yEdgeTemp(vectorIndex(obj.fineGrid,I2,J,K,'yedge'))'], ...
                    n*4,1)];
                entries = [ entries ; reshape( ...
                    [-ones(1,n) ; ones(1,n) ; ones(1,n) ; -ones(1,n)], ...
                    n*4,1)];
                %   some additional x- and y-edges need to be included, to allow
                %   for subdivided edges bounding a larger cell
                %  THIS DOES NOT ARISE FOR THE SPECIAL CASE WHERE ALL CELLS
                %  IN A LAYER HAVE THE SAME COARSNESS: for now I am giving
                %  up on generality!!!!
             ncol = length(obj.indE);
             nrow = length(obj.indF);
             obj.T = sparse(rowIndicies,colIndicies,entries,nrow,ncol);
        end
    end
end