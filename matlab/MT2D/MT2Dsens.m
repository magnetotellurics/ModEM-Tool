classdef MT2Dsens < Sensitivity
    % class for sensitivity matrix manipulations for 2D MT.
    % This version is based on full calculation of the matrix by call to
    % ModEM fortran code.  Idea is to use fixed methods names (i.e., make
    % this an instance of an abstract class)
    
    properties
        Sens   %   cell array of model parameters, one for each data
    end
    
    methods
        
        function obj = MT2Dsens(d,sigma0)
        %   class constructor.  With no arguments just makes an empty
        %    sensitivity object; with data vector input create empty cell
        %    array, with both d, sigma replicate function of makeSens.m,
        %    viz:
        %   Write sigma0, d to files to be read by Fortran program Mod2DMT
        %   calls program with option to  compute sensitivity matrix, reads result
        %    and returns result in J
        %
        %  Usage: J = makeSens(sigma0,d);
        %
        %  Inputs:	sigma0 = conductivity parameter (structure)
        %		d = data vector (cell array of structures)
            obj.dataVecClass = 'MT2DZ';
            obj.modelVecClass = 'MT2DmodelParam';
            obj.J_stored = true;
            
            if nargin > 0            
                if  ~isa(d,'MT2DZ')
                    error('data vector  of incorrect class')
                end
                obj.dataVecClass = class(d);
                obj.d = d;
                obj.N = lengthDat(d);
                obj.Sens = cell(obj.N,1);
            end
                
            if nargin == 2 
                obj.modelVecClass = class(sigma0);
                updateSens(obj,sigma0);
            end
        end
        %******************************************************************
        function updateSens(obj,m)
        %    upate sensitivity object for model parameter m
        %
        %    Usage: updateSens(obj,m)
        %
        %    Assumes object has been created and data vector structure has
        %          been set
        %    In this  implementation of the sensitivity matrix object 
        %     (MT2Dsens, full calculation of J) ModEM code is called to
        %      update full Jacobian
        
            if ~isa(m,'MT2DmodelParam') 
                error('model parameter of incorrect class')
            end

            %    make scratch directory if it doesn't exist ...
            if exist('scratch','dir')==0
                mkdir('scratch');
            end

            %   write out conductivity parameter in scratch directory
            m_File = 'scratch/Input.cpr';
            writeVec(m,m_File);

            %   write out data vector template in scratch directory
            d0_File = 'scratch/Input.imp';
            writeVec(obj.d,d0_File);
            %   file name for computed sensitivity matrix, to be output by Mod2DMT
            J_File = 'scratch/Out.sns';

            %   run Mod2DMT
            Test2D('COMPUTE_J',m_File,d0_File,J_File)

            obj.grid = m.grid;
            obj.paramType = m.paramType;
            obj.M = obj.grid.Ny*(obj.grid.Nz-obj.grid.Nza);
            obj.normalized = false;
            if obj.d.normalized
                obj.d = UnNormalize(obj.d);
            end
            
            %  read in sensitivity matrix
            readSens(obj,J_File);
            %obj.Sens = readCondMTX_2D(J_File);
     
        end
        %******************************************************************
        function readSens(obj,cfile)
            %
            % Usage:  [readSens(obj,cfile)
            %
            %  Reads in sequence of data sensitivities, return as cell
            %   array of structures ... one for each
            %   real observation.    These are not stored as modelParameter
            %   objects!
            
            fid = fopen(cfile,'r','n');
            
            fread(fid,1,'long');
            fread(fid,80,'char');
            fread(fid,2,'long');
            nSigma = fread(fid,1,'long');
            fread(fid,1,'long');
            
            for k = 1:nSigma
                fread(fid,1,'long');
                s = fread(fid,80,'char');
                s = deblank(char(s'));
                fread(fid,2,'long');
                nm = fread(fid,2,'long');
                ny = nm(1);
                nz = nm(2);
                fread(fid,2,'long');
                v = fread(fid,[ny,nz],'float64');
                fread(fid,2,'long');
                AirCond = fread(fid,1,'float64');
                fread(fid,1,'long');           
                obj.Sens{k} = struct('v',v,'paramType',s,'AirCond',AirCond);
            end
            fclose(fid);
        end
        %******************************************************************
        function dTilde = ExtractNormVec(obj,dIn)
        %  wrapper for ExtractVec (for MT2DZ objects just normalizes and
        %  extracts;  input "obj" is used to cause this instance of 
        %   ExtractNormVec to be called;
        %    same method for projected sensitivity (ProjSensMTX objects)
        %      also projects onto data subspace
           dTilde = ExtractVec(Normalize(dIn,1));
        end

        %******************************************************************
        function b = SetNormVec(obj,bTilde)
        %  wrapper for SetVec (for MT2DZ objects maps back to
        %     un-normalized data space--note that in the context of dataspace
        %       solution vectors this requires multiplication by CdInv
        %  Same method for projected sensitivity also completes the 
        %          back-projection)
           b = SetVec(obj.d,bTilde);
           if ~obj.normalized
               %    map to physical space to match sensitivity ... yes
               %    multiplying by CdInv is correct ...
               b =  MultCdInv(b,obj.d);
               %    need to reset the normalized attribute!
               b.normalized = false;
           else
               b.normalized = true;
           end
        end
        %******************************************************************
        function [dOut] = J_times_m(obj,m)
        %  Usage : [dOut] = J_times_m(J,m);
        %    multiply m by sensitivity matrix, put result in data space
        %    object dOut
 
%            if obj.normalized
%                error('')
%            end
            dOut = obj.d;
            kk = 0;
            for k = 1:dOut.NTX
                nSite = length(dOut.d{k}.siteLoc);
                for j = 1:nSite
                    kk = kk + 1;
                    Zr = sum(sum(obj.Sens{kk}.v .* m.v));
                    if dOut.d{k}.Cmplx
                        kk = kk + 1;
                        Zi = sum(sum(obj.Sens{kk}.v .* m.v));
                        dOut.d{k}.Z(j)  = Zr+1i*Zi;
                    else
                        dOut.d{k}.Z(j)  = Zr;
                    end
                end
            end
        end
        %******************************************************************
        function [mOut] = JT_times_d(obj,d)
        %  Usage : [mOut] = JT_times_d(J,d);
        % multiply by transpose of sensitivity matrix (stored as
        %     cell array of data/model space structures)
         
            if obj.normalized && ~d.normalized
                error('normalization of J and d inconsistent')
            end
            nTx = d.NTX;
            nJ = obj.d.NTX;
            if nJ ~= nTx  || length(obj.Sens) ~= lengthDat(d)
                error('Error: data vector and sensitivity size not consistent')
            end
            kk = 0;
            mTemp = obj.Sens{1};
            mTemp.v = zeros(size(mTemp.v));
            for k = 1:nTx
                nSite = length(d.d{k}.siteLoc);
                for j = 1:nSite
                    kk = kk + 1;
                    mTemp.v = mTemp.v + real(d.d{k}.Z(j))*obj.Sens{kk}.v;
                    if obj.d.d{k}.Cmplx
                        kk = kk + 1;
                        mTemp.v = mTemp.v + imag(d.d{k}.Z(j))*obj.Sens{kk}.v;
                    end
                end
            end
            mOut = MT2DmodelParam(obj.grid,obj.paramType);
            mOut = setModelParam(mOut,mTemp.v,mTemp.paramType,mTemp.AirCond);          
        end
        %******************************************************************
        function m =  rowJ(obj,k)
        %   extract row k from representer
            m = obj.Sens{k};
        end
        %******************************************************************
        function obj2 = copy(obj1)
        %   makes a copy of actual sensitivity object (not handle)
        
            obj2 = MT2Dsens;
            obj2.N = obj1.N;
            obj2.M = obj1.M;
            obj2.paramType = obj1.paramType;
            obj2.d = obj1.d;
            obj2.Sens = obj1.Sens;
            obj2.grid = obj1.grid;
        end
        %******************************************************************
        function Normalize(obj,CmHalf)
        %   Normalizes Jacobian structure, overwritting input Jacobian by
        %   J_tilde = C_m^{1/2} J C_d^{-1/2} and d_tilde = C_d{-1/2} d
            
            if obj.normalized
                error('Jacobian is already normalized')
            end
                       
            CdInv = InvDataErr(obj.d);
            for k = 1:obj.N
                obj.Sens{k} = CovMult(CmHalf,obj.Sens{k},1);
                obj.Sens{k}.v = obj.Sens{k}.v*CdInv(k);
            end
            %   also normalize stored data vector object ...
            obj.d = Normalize(obj.d);
            obj.normalized = true;
        end
        %******************************************************************    
        function [g] = JT_times_d_MTX(obj,d,SaveCmplx)
        %  variant on JT_times_d ... does not sum over
        %   frequencies, and optionally returns "imaginary part",
        %   for transmitters in input data vector d which are complex
        %   output is a cell array of size NTX or (NTX,2) for real and 
        %   complex cases, respectively
            
            if nargin < 3
                SaveCmplx = false;
            end    
            if obj.d.NTX ~= d.NTX
                error('Error: data vector and sensitivity size not consistent')
            end 
            if SaveCmplx
                g = cell(d.NTX,2);
            else
                g = cell(d.NTX,1);
            end
            kk = 0;
            modTemplate = MT2DmodelParam(obj.grid,obj.paramType);
            for k = 1:d.NTX
                mReal = zeros(size(obj.Sens{1}.v));
                nSite = length(d.d{k}.siteLoc);
                if SaveCmplx
                    mImag = mtemp;
                end
                for j = 1:nSite
                    kk = kk + 1;
                    mReal = mReal + real(d.d{k}.Z(j))*obj.Sens{kk}.v;
                    if SaveCmplx
                        mImag = mImag + imag(d.d{k}.Z(j))*obj.Sens{kk}.v;
                    end
                    
                    if obj.d.d{k}.Cmplx
                        kk = kk + 1;
                        mReal = mReal + imag(d.d{k}.Z(j))*obj.Sens{kk}.v;
                        if SaveCmplx
                            mImag = mImag - real(d.d{k}.Z(j))*obj.Sens{kk}.v;
                        end
                    end
                end
                g{k,1} = setModelParam(modTemplate,...
                    mReal,modTemplate.paramType,modTemplate.AirCond);
                if SaveCmplx
                    g{k,2} = setModelParam(modTemplate,...
                        mImag,modTemplate.paramType,modTemplate.AirCond);
                end              
            end
        end
    end    %  methods ....
        
end    %   classdef
