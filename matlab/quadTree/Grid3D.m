classdef Grid3D < handle
    % base class for storing 2D MT data; follows structure used in
    % original non-object-oriented "dataSpace" routines; idea is to
    % quickly develop a data-space object that can later be replaced by
    % something more logical
    
    properties
        Nx
        Ny    %     number of grid cells in y direction
        Nz    %     number of grid cells in z direction  (total, including air)
        Nza   %     number of air layers
        Dx    %     cell dimensions: x-direction
        Dy    %     cell dimensions: y-direction
        Dz    %     cell dimensions: z-direction
    end
    
    methods
        %*******************************************************************
        function obj = Grid3D(Dx,Dy,Dz,Nza)
            %   class constructor ... simple
            if nargin == 4
                obj.Nx = length(Dx);
                obj.Ny = length(Dy);
                obj.Nz = length(Dz);
                obj.Nza = Nza;
                obj.Dx = Dx;
                obj.Dy = Dy;
                obj.Dz = Dz;
            end
        end
        %*******************************************************************
        function gridOut = refineGridQT(obj,Kmax)
            %   refine grid (horizontal directions only) by a factor of
            %   2^Kmax;
            %   
            %   Usage: fineGrid = refineGridQT(coarseGrid,Kmax)
            %
            Nx = obj.Nx*2^Kmax;
            Dx = reshape(ones(2^Kmax,1)*obj.Dx'/2^Kmax,Nx,1);
            Ny = obj.Ny*2^Kmax;
            Dy = reshape(ones(2^Kmax,1)*obj.Dy'/2^Kmax,Ny,1);
            gridOut = Grid3D(Dx,Dy,obj.Dz,obj.Nza);
        end
        %*******************************************************************
        function [delX,delY,delZ] = dualLengths(obj)
            %   compute dual grid length elements
            delX = ([obj.Dx; 0]+[0; obj.Dx])/2;
            delY = ([obj.Dy; 0]+[0; obj.Dy])/2;
            delZ = ([obj.Dz; 0]+[0 ;obj.Dz])/2;
        end
        %*******************************************************************
        function [nx,ny,nz] = setLimits(obj,type)
            %   set arraay limits for components of types
            %                     cell, node,
            %                     xedge, yede, zedge
            %                     xface, yface, zface
             %   should probably add range checks
            switch lower(type)
                case 'cell'
                    nx = obj.Nx;
                    ny = obj.Ny;
                    nz = obj.Nz;
                case 'node'
                    nx = obj.Nx+1;
                    ny = obj.Ny+1;
                    nz = obj.Nz+1;
                case 'xedge'
                    nx = obj.Nx;
                    ny = obj.Ny+1;
                    nz = obj.Nz+1;
                case 'yedge'
                    nx = obj.Nx+1;
                    ny = obj.Ny;
                    nz = obj.Nz+1;
                case 'zedge'
                    nx = obj.Nx+1;
                    ny = obj.Ny+1;
                    nz = obj.Nz;
                case 'xface'
                    nx = obj.Nx+1;
                    ny = obj.Ny;
                    nz = obj.Nz;
                case 'yface'
                    nx = obj.Nx;
                    ny = obj.Ny+1;
                    nz = obj.Nz;
                case 'zface'
                    nx = obj.Nx;
                    ny = obj.Ny;
                    nz = obj.Nz+1;
            end
        end
        %*******************************************************************
        function [I,J,K] = gridIndex(obj,index,type)
            %  cells, nodes, edges, faces are enumerated as column vector
            %   elements in the standard way, consistent with 
            %    matlab reshape(X,Nx,Ny,Nz), where X is an array of
            %  dimension X(Nx,Ny,Nz).   Given array of integer indices of 
            %   vector elements this function, computes the
            %  corresponding i,j,k for an array of given type
            %   allowabel types : cell, node,
            %                     xedge, yede, zedge
            %                     xface, yface, zface
            
            [nx,ny,nz] = setLimits(obj,type);
            I = mod(index,nx);
            I(I==0) = nx;
            J  = mod(ceil(index/nx),ny);
            J(J==0) = ny;
            K = ceil(index/(nx*ny));
        end
        %*******************************************************************
        function [index] = vectorIndex(obj,I,J,K,type)
            %  cells, nodes, edges, faces are enumerated as column vector
            %   elements in the standard way, consistent with 
            %    matlab reshape(X,Nx,Ny,Nz), where X is an array of
            %  dimension X(Nx,Ny,Nz).   Given array of integer indices of 
            %   vector elements this function, computes the
            %  corresponding i,j,k for an array of given type
            %   allowabel types : cell, node,
            %                     xedge, yede, zedge
            %                     xface, yface, zface
            
            [nx,ny,nz] = setLimits(obj,type);
            index = (K-1)*nx*ny+(J-1)*nx+I;
        end
    end     % methods
end    % classdef