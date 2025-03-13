classdef SensitivityIndirect < Sensitivity
    % class for sensitivity matrix manipulations , using indirect
    %    approach (i.e., not calculating the full sensitivity first).
    % This version is a variant on the version which uses a full calculation
    %     of the matrix, also using calls to
    % ModEM fortran code.  Idea is to use fixed methods names (i.e., make
    % this an instance of an abstract class)
    
    properties
        m0     %  (background) model parameter used to calculate sensitivity
        Cm     %   pointer to covariance object
    end
    
    methods
        
        function obj = SensitivityIndirect(d,m)
            %   class constructor
            %
            %  Inputs:	sigma = conductivity parameter (structure)
            %		d = data vector (cell array of structures)
            obj.dataVecClass = class(d);
            obj.d = d;
            obj.N = lengthDat(d);
            obj.m0 = m;
            obj.M = ModelParamLength(m);
            obj.modelVecClass = class(m);
            obj.J_stored = false;
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
            
            if ~isa(m,obj.modelVecClass)
                error('model parameter of incorrect class')
            end
            
            obj.m0 = m;
            obj.normalized = false;
            
        end
        %******************************************************************
        function dTilde = ExtractNormVec(obj,dIn) %#ok<MANU>
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
            %    make scratch directory if it doesn't exist ...
            if exist('scratch','dir')==0
                mkdir('scratch');
            end
            
            %   write scratch file for input conductivity parameter
            m0_File = 'scratch/m0.cpr';
            writeVec(obj.m0,m0_File);
            %   write scratch file for input data vector (template)
            d0_File = 'scratch/d0.imp';
            writeVec(obj.d,d0_File);
            %   write scratch file for input conductivity parameter perturbation
            dm_File = 'scratch/dm.cpr';
            if obj.normalized
                %   multiply input model parameter by C_m^{1/2}
                m = CovMult(obj.Cm,m,1);
            end
            writeVec(m,dm_File);
            %   file for output data, to be read after running program
            Jxdm_File = 'scratch/Jxdm.imp';
            
            %  execute program with -G option to multiply by sensitivity matrix
            if strcmp(obj.dataVecClass,'MT2DZ')
                Test2D('MULT_BY_J',m0_File,dm_File,d0_File,Jxdm_File)
            elseif strcmp(obj.dataVecClass,'MT3DZ')
                Test3D('MULT_BY_J',m0_File,dm_File,d0_File,Jxdm_File)
            end
            [dOut] = readVec(obj.d,Jxdm_File);
            
            if obj.normalized
                dOut = MultCdInv(dOut,obj.d);
            end
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
            
            %    make scratch directory if it doesn't exist ...
            if exist('scratch','dir')==0
                mkdir('scratch');
            end
            m0_File = 'scratch/m0.cpr';
            writeVec(obj.m0,m0_File);
            %   write scratch file for input data vector
            d_File = 'scratch/d.imp';
            writeVec(d,d_File);
            %  file name for output sensitivity (model parameter)
            JTxd_File = 'scratch/JTxd.cpr';
            
            if strcmp(obj.dataVecClass,'MT2DZ')
                Test2D('MULT_BY_J_T',m0_File,d_File,JTxd_File)
            elseif strcmp(obj.dataVecClass,'MT3DZ')
                Test3D('MULT_BY_J_T',m0_File,d_File,JTxd_File)
            end
            mOut = readVec(obj.m0,JTxd_File);
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
            %   here m is returned as a model parameter object ... note
            %   that not all versions of rowJ are not returning m as an
            %   actual model parameter object!   Only m.v is actually used,
            %   so we are "getting away" with this!
        end
        %******************************************************************
        function obj2 = copy(obj1)
            %   makes a copy of actual sensitivity object (not handle)
            
            obj2 = SensitivityIndirect(obj1.d,obj1.m0);
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
            
            if nargin < 3
                SaveCmplx = false;
            end
            if SaveCmplx
                error('Cant comppute imaginary parts of J^Td at present')
            end
            if obj.normalized
                if d.normalized==0
                    error('normalization of J and d inconsistent')
                elseif d.normalized==1
                    d = MultCdInv(d,obj.d);
                end
            end
            
            %    make scratch directory if it doesn't exist ...
            if exist('scratch','dir')==0
                mkdir('scratch');
            end
            m0_File = 'scratch/m0.cpr';
            writeVec(obj.m0,m0_File);
            %   write scratch file for input data vector
            d_File = 'scratch/d.imp';
            writeVec(d,d_File);
            %  file name for output sensitivity (model parameter)
            JTxd_File = 'scratch/JTxd.mtx';
            
            if strcmp(obj.dataVecClass,'MT3DZ')
                Test3D('MULT_BY_J_T_Mtx',m0_File,d_File,JTxd_File)
                g = readCondMTX_3D(JTxd_File);
            else
                error('Indirect version of JmultT_MTX only avaialble for 3D')
            end
            nTx = length(g);
            for iTx = 1:nTx
                if obj.normalized
                    g{iTx} = CovMult(obj.Cm,g{iTx},1);
                end
            end
        end
    end    %  methods ....
    methods (Static)
        function [status] = Run2D(todo,arg1,arg2,arg3,arg4,arg5,arg6)
            
            if nargin == 3
                args = [arg1 ' ' arg2];
            elseif nargin == 4
                args = [arg1 ' ' arg2 ' ' arg3];
            elseif nargin == 5
                args = [arg1 ' ' arg2 ' ' arg3 ' ' arg4];
            elseif nargin == 6
                args = [arg1 ' ' arg2 ' ' arg3 ' ' arg4 ' ' arg5];
            elseif nargin == 7
                args = [arg1 ' ' arg2 ' ' arg3 ' ' arg4 ' ' arg5 ' ' arg6];
            else
                error('Wrong number of arguments to Run2D');
            end
            
            switch todo
                case 'READ_WRITE'
                    disp(['Mod2DMT -R ' args]);
                    status = system(['Mod2DMT -R ' args]);
                case 'FORWARD'
                    disp(['Mod2DMT -F ' args]);
                    status = system(['Mod2DMT -F ' args]);
                case 'COMPUTE_J'
                    disp(['Mod2DMT -J ' args]);
                    status = system(['Mod2DMT -J ' args]);
                case 'MULT_BY_J'
                    disp(['Mod2DMT -M ' args]);
                    status = system(['Mod2DMT -M ' args]);
                case 'MULT_BY_J_T'
                    disp(['Mod2DMT -T ' args]);
                    status = system(['Mod2DMT -T ' args]);
                case 'MULT_BY_J_T_MTX'
                    disp(['Mod2DMT -x ' args]);
                    status = system(['Mod2DMT -x ' args]);
                case 'INVERSE_NLCG'
                    disp(['Mod2DMT -I NLCG ' args]);
                    status = system(['Mod2DMT -I NLCG ' args]);
                otherwise
                    error('This job is not implemented');
            end
        end
        %******************************************************************
        function [status] = Run3D(todo,arg1,arg2,arg3,arg4,arg5,arg6)
            
            if nargin == 3
                args = [arg1 ' ' arg2];
            elseif nargin == 4
                args = [arg1 ' ' arg2 ' ' arg3];
            elseif nargin == 5
                args = [arg1 ' ' arg2 ' ' arg3 ' ' arg4];
            elseif nargin == 6
                args = [arg1 ' ' arg2 ' ' arg3 ' ' arg4 ' ' arg5];
            elseif nargin == 7
                args = [arg1 ' ' arg2 ' ' arg3 ' ' arg4 ' ' arg5 ' ' arg6];
            else
                error('Wrong number of arguments to Test3D');
            end
            
            switch todo
                case 'READ_WRITE'
                    disp(['Mod3DMT -R ' args]);
                    status = system(['Mod3DMT -R ' args]);
                case 'FORWARD'
                    disp(['Mod3DMT -F ' args]);
                    status = system(['Mod3DMT -F ' args]);
                case 'COMPUTE_J'
                    disp(['Mod3DMT -J ' args]);
                    status = system(['Mod3DMT -J ' args]);
                case 'MULT_BY_J'
                    disp(['Mod3DMT -M ' args]);
                    status = system(['Mod3DMT -M ' args]);
                case 'MULT_BY_J_T'
                    disp(['Mod3DMT -T ' args]);
                    status = system(['Mod3DMT -T ' args]);
                case 'MULT_BY_J_T_MTX'
                    disp(['Mod3DMT -X ' args]);
                    status = system(['mpirunMod2DMT -x ' args]);
                case 'INVERSE_NLCG'
                    disp(['Mod3DMT -I NLCG ' args]);
                    status = system(['Mod3DMT -I NLCG ' args]);
                otherwise
                    error('This job is not implemented');
            end
        end
    end
    
end    %   classdef
