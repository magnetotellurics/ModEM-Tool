classdef ProjSensMTX < Sensitivity
    %  Class for projectd normalized sensitivity matrix manipulations,
    %     explicitly allowing for multiple transmitters
    %  Projected data vectors (which J^T multiplies) are just standard real
    %    vectors.
    %
    %   Gary D. Egbert, 2010
    %   College of Oceanic and Atmospheric Sciences
    %   Oregon State University
    
    properties
        Sens   %    cell array of model parameters arrays, one for each transmitter
        NTx    %    number of transmitters
        Kmax   %    maximum number of sensitivities for one transmitter
        K      %     number of sensitivity contrasts (same for each transmitter)
        U      %    cell array of K (normalized) data vectors ... in 
               %     principal these are orthonormal, but they might not
               %     be!
        B      %    array of size(2,N+1) used to store bidagonal matrix from
               %         Lancozs process
        truncErr    %    estimate of truncation error 
        jk     %    index array to map data vector components to transmitters
               %     only used when MTX (mult-transmitter) option is true
        MTX = false    %  true to save sensitivities for individul
               %     transmitters
    end
    
    methods     
        function obj = ProjSensMTX(NTx,Kmax,saveMTX)
            if nargin == 3
                obj.MTX =  saveMTX;
            end
            %   simple class constructor.
            obj.Kmax = Kmax;
            obj.NTx = NTx;
            obj.U = cell(Kmax,1);
            obj.B = zeros(2,Kmax);
            obj.truncErr = zeros(Kmax,1);
            obj.normalized = true;
            if obj.MTX
                obj.Sens = cell(NTx,Kmax);
                obj.jk = zeros(2,Kmax*NTx);
            else
                obj.Sens = cell(1,Kmax);
            end
            obj.N = 0;
            obj.K = 0;
        end
        %******************************************************************
        function updateSens(obj,U1,Sens1,beta_alpha)
        %    adds single data vector and sensitivities (possibly for each
        %    transmitter separately) corresponding to this linear
        %    combination of data to  ProjSensMTX object
        %    Now assuming there is exactly one component f
        %    transmitter
        
            if obj.N == 0
                %   some things to initialize still ...
                obj.dataVecClass = class(U1);
                if obj.MTX
                    obj.modelVecClass = class(Sens1{1});
                    obj.M = ModelParamLength(Sens1{1});
                else
                    obj.modelVecClass = class(Sens1);
                    obj.M = ModelParamLength(Sens1);
                end
            end
            
            if obj.MTX
                n = obj.K*obj.NTx;
                obj.K = obj.K+1;
                for j = 1:obj.NTx
                    n = n + 1;
                    obj.Sens{j,obj.K} = Sens1{j};
                    obj.jk(1,n)  = j;
                    obj.jk(2,n) = obj.K;
                end
                obj.N = n;
            else 
                obj.K = obj.K+1;
                obj.Sens{obj.K} = Sens1;
                obj.N = obj.K;
            end
            obj.U{obj.K} = U1;
            if(nargin>3)
                obj.B(1,obj.K) = beta_alpha(1);
                obj.B(2,obj.K) = beta_alpha(2);
                obj.truncErr(obj.K) = beta_alpha(3);
            end
        end
        %******************************************************************
        function [dOut] = J_times_m(obj,m)
        %  Usage : [dOut] = J_times_m(J,m);
        %    multiply m by projected normalized sensitivity matrix,
        %    result is a standard real matlab vector
            dOut = zeros(obj.N,1);
            if obj.MTX
                for n = 1:obj.N
                    dOut(n) = dot(obj.Sens(obj.jk(1,n),obj.jk(2,n)),m);
                end
            else
                for n = 1:obj.N
                    dOut(n) = dot(obj.Sens(n),m);
                end
            end
        end
        %******************************************************************
        function [mOut] = JT_times_d(obj,d)
        %  Usage : [mOut] = JT_times_d(J,d);
        %  multiply by transpose of projected normalized sensitivity matrix 
        %    (stored as cell array of model space structures)
        %  d is a standard real matlab vector ... implicitly assumed to be
        %     normalized!
            Nd = length(d);
            if Nd ~= obj.N  
                error('Error: data vector and sensitivity size not consistent')
            end
            n = 0;
            mOut = ZeroVec(obj.Sens{1});
            for k = 1:obj.K
                if obj.MTX
                    for j = 1:obj.NTx
                        n = n + 1;
                        mOut = mOut + d(n)*...
                            obj.Sens{obj.jk(1,n),obj.jk(2,n)};
                    end
                else
                    mOut = mOut + d(k)*obj.Sens{k};
                end
            end    
        end
        %******************************************************************
        function m =  rowJ(obj,n)
        %   extract row n from representer
            if obj.MTX
               m = obj.Sens{obj.jk(1,n),obj.jk(2,n)};
            else
               m = obj.Sens{n};
            end
        end
        %******************************************************************
        function dTilde = ExtractNormVec(obj,dIn)
        %  for projected sensitivity normalizes dIn, and then computes
        %  linear combinations stored in obj.U
           dIn = Normalize(dIn,1);
           dTilde = zeros(obj.N,1);
           for k = 1:obj.K
              if obj.MTX
                  n1 = (k-1)*dIn.NTX+1;
                  n2 = k*dIn.NTX;
                  dTilde(n1:n2) = dotMTX(obj.U{k},dIn);
              else 
                  dTilde(k) = dot(obj.U{k},dIn);
              end
           end
        end
        %******************************************************************
        function b = SetNormVec(obj,bTilde)
        %   maps back to the original data space; this simple version
        %    assumes that obj.U is exactly orthonogonal.   Need to
        %     orthogonalize *before* doing any RepMat operations
        %   
        %   After converting back to a data space object, the
        %   result is multiplied by CdInv to map to coefficients for
        %    unnormalized problem
        
            b = zeroVec(obj.d);
            b.normalized = obj.normalized;
            n =  0;
            for k = 1:obj.K
               if obj.MTX
                   for j = 1:obj.NTx
                      n = n+1;
                      b = add1Tx(b,j,bTilde(n),obj.U{k});
                   end
               else
                   b = b + bTilde(k)*obj.U{k};
               end
            end
            if ~obj.normalized
                b =  MultCdInv(b);
                b.normalized = false;
            else
                b.normalized = true;
            end
        end
        %******************************************************************
        %
        function [m,dRes,relErr] = LanczosStep(obj,J,dTilde,nu)
        %    Given an input ProjSensMTX object (obj), the orignal
        %     normalized sensitivity matrix object (J), and the (RHS vector
        %      in the normal equations, make the next step in the
        %    "generalized Lanczos process"   ... what to call this routine?
        
            if ~dTilde.normalized
                dTilde = Normalize(dTilde);
            end
            if obj.K > 0
                %   start by solving projected problem
                R = RepMat(obj);
                R.formR
                R.factorR
                b = solveR(R,dTilde,nu);
                %   R.btilde is the solution in the projected (and
                %   normalized) space ... output b is in original 
                %   data space   ... normalized if R.J is (i.e., if obj is)
                m = JT_times_d(obj,R.bTilde);
                %     now m is in the transformed model space ...  and
                %    note that the J_times_M used below is NOT the
                %          method in this class, but is for the original 
                %     (unprojected) sensitivity matrix class
                %     J should  be normalized here ... so dFit will be also
                dFit = J_times_m(J,m);
            else
                %   start from zero
                dFit = zeroVec(dTilde);
                m = [];
            end
            dRes = dTilde-(dFit+nu*b);        
            relErr = sqrt(dot(dRes,dRes)/(dot(dTilde,dTilde)));       
        end
        %******************************************************************
        function [m,ITER] = Lanczos(obj,J,dHat,CmHalf,ITER,nu)
        %    this function does a modified Lanczos bi-diagonalization,
        %    filling out the projected sensitivity matrix obj   ... this
        %    "ProjSensMTX" object should be created before calling.
        %    Optionally, some data contrasts/sensitivities could already be
        %    stored in obj before calling.  And if obj.MTX = true the
        %    multi-transmitter variant on the Lanczos scheme is used.
            
            %   make sure input Jacobian is normalized
            if ~J.normalized
                Normalize(J,CmHalf);
            end
            if ~dHat.normalized
                dTilde = Normalize(dHat,1);
            else
                dTilde = dHat;
            end
            if obj.N > 0
                %  some sensitivities already calculated ... solve linear
                %  problem in projected space to compute current data space
                %  residual
                [m,dRes,ITER.relErr(1)] = LanczosStep(obj,J,dTilde,nu);
            else
                %   just start from input (ater normalization)
                dRes = dTilde;
                ITER.relErr(1) = 1;
            end
            
            DONE = false;
            k=1;
            while ~DONE
                k=k+1;
                U1 = UnitVec(dRes,obj.MTX);
                if obj.MTX
                    %  here m will be an array of NTX model parameters
                    m = JT_times_d_MTX(J,U1);
                else
                    m = JT_times_d(J,U1);
                end
                updateSens(obj,U1,m);
                [m,dRes,ITER.relErr(k)] = LanczosStep(obj,J,dTilde,nu);
                ITER.niter = k;
                DONE = checkConvergence(ITER);
            end
        end
        %******************************************************************
        function [dVecs] = significant_dVec_MTX(obj,lambdaMin)
        %  given ProjSensMTX object, find significant data vectors (based
        %  on eigenvector decomposition of RepMat constructed from object)
        %  Usage: [dVecs] = significant_dVec_MTX(obj,lambdaMin)
            if ~obj.MTX
                error('significant_dvec_MT can only be used for MTX case')
            end
            R = RepMat(obj);
            R.formR
            R.factorR
            nUse = find((R.lambda < lambdaMin),1)-1;
            dVecs = cell(nUse,1);
            for iUse = 1:nUse
                dVecs{iUse} = zeroVec(obj.U{1});
                for k = 1:obj.K
                    i1 = (k-1)*obj.NTx;
                    i2 = k*obj.NTx;
                    c = conj(R.U(i1:i2,iUse));
                    dVecs{iUse} = addMTX(dVecs{iUse},k,c,obj.U{k});
                end        
                
                
            end        
        end
        %******************************************************************
        function FillProjSensMTX(obj,dVecs,J)
        %   given a set of data space vectors and a sensitivity matrix
        %   object J, fill in ProjSensMTX object (obj)     
            for k = 1:length(dVecs)
                 %  here m will be an array of NTX model parameters
                 m = JT_times_d_MTX(J,dVecs{k});
                 updateSens(obj,dVecs{k},m)
            end
        end
        %******************************************************************
        function Orthogonalize(obj)
        %   Simple explicit orthonormalization of obj.U, corresponding
        %   transformation of obj.Sens ... 
        
            if obj.MTX
                 a = 1./sqrt(dotMTX(obj.U{1},obj.U{1}));
                 for l = 1:obj.U{1}.NTX
                     obj.U{1} =  mult1Tx(obj.U{1},l,a(l));
                     obj.Sens{l,1} = a(l)*obj.Sens{l,1};
                 end
                for k = 2:obj.K
                    for j = 1:k-1
                        gamma = dotMTX(obj.U{j},obj.U{k});
                        for l = 1:obj.U{1}.NTX
                            obj.U{k} = add1Tx(obj.U{k},l,-gamma(l),obj.U{j});
                            obj.Sens{l,k} = obj.Sens{l,k}-gamma(l)*obj.Sens{l,j};
                        end
                    end
                    a = 1./sqrt(dotMTX(obj.U{k},obj.U{k}));
                    for l = 1:obj.U{1}.NTX
                        obj.U{k} =  mult1Tx(obj.U{k},l,a(l));
                        obj.Sens{l,k} = a(l)*obj.Sens{l,k};
                    end
                end
            else
                a = 1./sqrt(dot(obj.U{1},obj.U{1}));
                obj.U{1} = a*obj.U{1};
                obj.Sens{1} = a*obj.Sens{1};
                for k = 2:obj.K
                    for j = 1:k-1
                        gamma = dot(obj.U{j},obj.U{k});
                        obj.U{k} = obj.U{k}-gamma*obj.U{j};
                        obj.Sens{k} = obj.Sens{k}-gamma*obj.Sens{j};
                    end
                    a = 1./sqrt(dot(obj.U{k},obj.U{k}));
                    obj.U{k} = a*obj.U{k};
                    obj.Sens{k} = a*obj.Sens{k};
                end
            end
        end
        %******************************************************************
        function [UU] = UprimeU(obj)
        %   compute U'U from data objects ... for checking orthonormality
            UU = zeros(obj.K);
            for j = 1:obj.K
                for k = 1:obj.K
                    UU(j,k) = dot(obj.U{j},obj.U{k});
                end
            end
        end
        %******************************************************************
        function [relErr] = JprojRMS(Jproj,J,dTilde,nu)
            
            %   given J, dTilde used to derive ProjSensMTX object
            %    compute RMS error in solution of linear normal equations
            %    for each step k= 1, ... , K
            
            if ~dTilde.normalized
                dTilde = Normalize(dTilde);
            end
             
            K = Jproj.K;
            relErr = zeros(K,1);
            
            for k = 1:K
                
                Jproj.K = k;
                if Jproj.MTX
                    Jproj.N = Jproj.NTx*k;
                else
                    Jproj.N = k;
                end
                %   start by solving projected problem
                R = RepMat(Jproj);
                R.formR
                R.factorR
                b = solveR(R,dTilde,nu);
                %   R.btilde is the solution in the projected (and
                %   normalized) space ... output b is in original
                %   data space   ... normalized if R.J is (i.e., if obj is)
                m = JT_times_d(Jproj,R.bTilde);
                %     now m is in the transformed model space ...  and
                %    note that the J_times_M used below is NOT the
                %          method in this class, but is for the original
                %     (unprojected) sensitivity matrix class
                %     J should  be normalized here ... so dFit will be also
                dFit = J_times_m(J,m);
                
                dRes = dTilde-(dFit+nu*b);
                relErr(k) = sqrt(dot(dRes,dRes)/(dot(dTilde,dTilde)));
            end 
        end
        %******************************************************************
        function [obj] = collapseMTX(objMTX)
            %  create Jproj object with MTX = false from  an object with
            %    MTX = true
            
            saveMTX=false;
            obj = ProjSensMTX(objMTX.NTx,objMTX.Kmax,saveMTX);
            obj.N = objMTX.K;
            obj.M = objMTX.M;
            obj.K = objMTX.K;
            obj.U  = objMTX.U;
            obj.d = objMTX.d;
            obj.B = objMTX.B;
            obj.truncErr = objMTX.truncErr;
            obj.dataVecClass = objMTX.dataVecClass;
            obj.modelVecClass = objMTX.modelVecClass;
            for k = 1:obj.K
                obj.Sens{k} = objMTX.Sens{1,k};
                for l = 2:obj.NTx
                    obj.Sens{k} = obj.Sens{k}+objMTX.Sens{l,k};
                end
            end
        end
        
        %    The following were created as dummy routines, needed because I
        %     made ProjSensMTX a sub-class of the abstract class
        %      Sensitivity, and these are all declared as abstract methods 
        %        for the superclass.
        %******************************************************************
        function obj2 = copy(obj1)
        %   makes a copy of actual sensitivity object (not handle)
        end

        %******************************************************************
        function Normalize(obj,CmHalf)
        %   Normalizes Jacobian structure, overwritting input Jacobian by
        %   J_tilde = C_m^{1/2} J C_d^{-1/2} and d_tilde = C_d{-1/2} d
        end

        %******************************************************************    
        function [g] = JT_times_d_MTX(obj,d,SaveCmplx)
        %  variant on JT_times_d ... does not sum over
        %   frequencies, and optionally returns "imaginary part",
        %   for transmitters in input data vector d which are complex
        %   output is a cell array of size NTX or (NTX,2) for real and 
        %   complex cases, respectively
        end
    end    %  methods ....
        
end    %   classdef
