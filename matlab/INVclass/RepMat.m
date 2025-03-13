classdef RepMat < handle
    % class for "representer matrix", essentialy
    %    CdInv * J * Cm^2 * J' * CdInv +lambda I
    %   coefficient matrix for data space normal equations
    
    properties
        J       %   sensitivity matrix handle
        nu=1    %   tradeoff parameter
        Cm      %   model covariance handle
        CdInv      %   data error covariance
        grid    %   This is just a handle of the grid object, making it
        storeR = false
        R       %    N x N storage for cross-products  ... only used if 
                %      storeR=true
        U       %    eigenvectors of R
        lambda  %    eigenvalues of R
        vecClass  %   class name for vectors that R should multiply
        dTilde  %  normalizeded, projected, extracted, etc. real N-vector
                %   that forms the rhs for the normal  equations
        bTilde  %  solution vector, again a simple real N-vector 
    end
    
    methods
        
        function obj = RepMat(J,Cm,nu)
        %   class  constructor
            if nargin >=1
                obj.J = J;
                obj.vecClass = J.dataVecClass;
                %   need to set data errors ... 
                obj.CdInv = InvDataErr(obj.J.d);
            end
            if nargin >=2
                obj.Cm = Cm;
            end
            if nargin == 3;
                obj.nu = nu;
            end
        end
        %******************************************************************
        function updateRepMat(obj,J)
        %   updates representer matrix, using new Jacobian
            obj.J = J;
        end
        %******************************************************************
        function x = mtimes(obj,b)
        %  compute matrix vector product x = (R+nuI)b
        
            if ~isa(b,obj.vecClass)
                error('matrix and vector of incompatible class for multiplication')
            end
            if ~b.normalized
                error('input vector x for product R*x should be normalized')
            end
            if obj.storeR
                temp = obj.R*ExtractVec(b);
                x = SetVec(b,temp)+obj.nu*b;
            else              
                if obj.J.normalized
                    y = JT_times_d(obj.J,b);
                    x = J_times_m(obj.J,y);
                    x = x+obj.nu*b;
                else
                    x = MultCdInv(b,obj.J.d);
                    y = JT_times_d(obj.J,x);
                    y = CovMult(obj.Cm,y,2);
                    x = J_times_m(obj.J,y);
                    x = MultCdInv(x,obj.J.d);
                    x = x+obj.nu*b;
                end
            end
        end
        %******************************************************************
        function formR(obj)
        %   makes full data-space cross-product matrix, normalized by 
        %    data error standard deviations
           obj.R = zeros(obj.J.N,obj.J.N);
           Jmtrx = zeros(obj.J.N,obj.J.M);
           
           for j = 1:obj.J.N
               mj = rowJ(obj.J,j);
               if ~obj.J.normalized
                   mj = CovMult(obj.Cm,mj,1);
               end
               %Jmtrx(j,:) = ExtractVec(mj);
               Jmtrx(j,:) = reshape(mj.v,1,obj.J.M);
           end
           if ~obj.J.normalized
               Jmtrx = diag(obj.CdInv)*Jmtrx;
           end
           
           obj.R = Jmtrx*Jmtrx';
           obj.storeR = true;
        end
        
        %******************************************************************        
        function factorR(obj)
        %  eigenvector decomposition of R
           [V,D] = eig(obj.R);
           [temp,I] = sort(real(diag(D)));
           obj.lambda = temp(end:-1:1);
           obj.lambda(obj.lambda<0) = 0;
           obj.U = V(:,I(end:-1:1));
        end
        %******************************************************************        
        function b = solveR(obj,dHat,nuIn)
        %  use factored representer matrix to solve (R + nu C_d)b = d
        %     the normalization state of the data-space vector b that is
        %     returned matches the normlization state of the input obj.J
           if nargin ==3
               obj.nu = nuIn;
           end
           obj.dTilde = ExtractNormVec(obj.J,dHat);
           obj.dTilde = obj.U'*obj.dTilde;
           obj.bTilde = obj.dTilde./(obj.lambda+obj.nu);
           obj.bTilde = obj.U*obj.bTilde;
           b = SetNormVec(obj.J,obj.bTilde);
           % SetUnNovmVec puts b back in physical data space if J is
           % unormzlized, in normalized data space if J is normalized
        end
        %******************************************************************
        function [nu] = SETnu(obj)
        %  Simple initialization of array of damping parameters
        %   to try in our crude search implementation
        %
        %  Usage:  [nu] = SETnu(lambda1)
            
            %   parameters for automatic generation of nu
            n_nu = 12;
            nu_1 = 0; nu_2 = 3.;
            lambda0 = mean(log10(obj.lambda(1:2)));
            
            nu1 = lambda0-nu_1;
            nu2 = nu1+nu_1-nu_2;
            dnu = (nu1-nu2)/n_nu;
            nu = nu2:dnu:nu1;
            nu = 10.^nu;
        end
        %******************************************************************
        function [S,V,U] = svdJ(R,lambda0)
        % makes truncated svd of (normalized) J using representer matrix R
        
           nTrunc = find(R.lambda<lambda0,1);
           U = R.U(:,1:nTrunc);
           Jmtrx = zeros(R.J.N,R.J.M);
           
           for j = 1:R.J.N
               mj = rowJ(R.J,j);
               if ~R.J.normalized
                   mj = CovMult(R.Cm,mj,1);
               end
               %Jmtrx(j,:) = ExtractVec(mj);
               Jmtrx(j,:) = reshape(mj.v,1,R.J.M);
           end
           if ~R.J.normalized
               Jmtrx = diag(R.CdInv)*Jmtrx;
           end
           
           V = U'*Jmtrx;
           S = sqrt(sum(V.*V,2));
           V = diag(1./S)*V;
           
        end
    end
end

