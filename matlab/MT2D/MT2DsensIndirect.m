classdef MT2DsensIndirect < Sensitivity
    % class for sensitivity matrix manipulations for 2D MT, using indirect
    %    approach (i.e., not calculating the full sensitivity first).
    % This version is a variant on the version which uses a full calculation 
    %     of the matrix, also using calls to
    % ModEM fortran code.  Idea is to use fixed methods names (i.e., make
    % this an instance of an abstract class)
    
    properties
        m0     %  (background) model parameter used to calculate sensitivity
        dataVecClass = 'MT2DZ'   %    data vector class for this
                                 %    implementation of sensitivity
        modelVecClass = 'MT2DmodelParam'   %  model parameter class
    end
    
    methods
        
        function obj = MT2DsensIndirect(d,m)
        %   class constructor.  With no arguments just makes an empty
        %    sensitivity object; otherwise initialize as much as possible
        %     given inputs available
        %
        %  Usage: J = makeSens(sigma,d);
        %
        %  Inputs:	sigma = conductivity parameter (structure)
        %		d = data vector (cell array of structures)
            
           if nargin > 0
                               
                if  ~isa(d,'MT2DZ')
                    error('data vector  of incorrect class')
                end
                
                obj.d = d;
                obj.N = lengthDat(d);
            end
                
            if nargin == 2            
                updateSens(obj,m);
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

           obj.m0 = m;
           obj.M = ModelParamLength(m);
     
        end
        %******************************************************************
        function [dOut] = J_times_m(obj,m)
        %  Usage : [dOut] = J_times_m(J,m);
        %    multiply m by sensitivity matrix, put result in data space
        %    object dOut
             %    make scratch directory if it doesn't exist ...
             if exist('scratch','dir')==0
                 mkdir('scratch');
             end
             %   write scratch file for input conductivity parameter
             m0_File = 'scratch/m0.cpr';
             writeVec(m0_File,obj.m0);
             %   write scratch file for input data vector (template)
             d0_File = 'scratch/d0.imp';
             writeVec(d0_File,obj.d);
             %   write scratch file for input conductivity parameter perturbation
             dm_File = 'scratch/dm.cpr';
             writeVec(dm_File,m);
             %   file for output data, to be read after running program
             Jxdm_File = 'scratch/Jxdm.imp';
             
             %  execute program with -G option to multiply by sensitivity matrix
             Test2D('MULT_BY_J',m0_File,dm_File,d0_File,Jxdm_File)
             [dOut] = readVec(obj.d,Jxdm_File);
             
        end
        %******************************************************************
        function [mOut] = JT_times_d(obj,d)
            %  Usage : [mOut] = JT_times_d(J,d);
            % multiply by transpose of sensitivity matrix (stored as
            %     cell array of data/model space structures)
            
            if obj.normalized 
                if d.normalized==0
                   error('normalization of J and d inconsistent')
                elseif d.normalized==1
                   d = MultCdInv(d,obj.d);
                end
            end
            nTx = d.NTX;
            nJ = obj.d.NTX;
            if nJ ~= nTx  || length(obj.Sens) ~= lengthDat(d)
                error('Error: data vector and sensitivity size not consistent')
            end
            %    make scratch directory if it doesn't exist ...
            if exist('scratch','dir')==0
                mkdir('scratch');
            end
            m0_File = 'scratch/m0.cpr';
            writeVec(m0_File,obj.m0);
            %   write scratch file for input data vector
            d_File = 'scratch/d.imp';
            writeVec(d_File,d);
            %  file name for output sensitivity (model parameter)
            JTxd_File = 'scratch/JTxd.cpr';
            
            Test2D('MULT_BY_J_T',m0_File,d_File,JTxd_File)
            mOut = readVec(obj.m0,JTxdFile);
            if obj.normalized
                mOut = CovMult(obj.Cm,mOut,1);
            end
        end
        %******************************************************************
        function m =  rowJ(obj,k)
        %   compute row of representer J^T e_k   where e_k is unit vector
        %   for component k of sensitvity
            eVec = zeros(obj.N,1);
            eVec(k) = 1;
            e = SetVec(obj.d,eVec);
            m = JT_times_d(obj,e);
        end
        %******************************************************************
        function obj2 = copy(obj1)
        %   makes a copy of actual sensitivity object (not handle)
        
            obj2 = MT2DsensIndirect;
            obj2.d = obj1.d;
            obj2.m0 = obj1.m0;
            obj2.grid = obj1.grid;
        end
        %******************************************************************
        function Normalize(obj,Cm)
        %   Normalizes Jacobian structure ... in indirect case just stores
        %   pointer to covariance and sets normalized = true
            
            obj.normalized = true;
            obj.Cm = Cm;
        end
        %******************************************************************    
        function [g] = JT_times_d_MTX(obj,d,SaveCmplx)
        %  variant on JT_times_d ... does not sum over
        %   frequencies, and optionally returns "imaginary part",
        %   for transmitters in input data vector d which are complex
        %   output is a cell array of size NTX or (NTX,2) for real and 
        %   complex cases, respectively
            
           error('JT_times_d_MTX not code for indirect sensitivity case')
        end
    end    %  methods ....
        
end    %   classdef
