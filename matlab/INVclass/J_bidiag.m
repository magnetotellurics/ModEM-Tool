function [Jproj] = J_bidiag(J,CmHalf,dHat,K,saveMTX)

%   Usage: [Jproj] = J_bidiag(J,d,K)
%   Computes K-step Lanczos bidagonalization of Cm^{1/2} J Cd^{-1/2}
%  implements bi-diag1 (Lower triangular)
%   decomposition of Paige and Saunders (1982).
%
%   if optional argument saveMTX=true save separate sensitivities for each
%   transmitter (defaults to false)
%
%   Note that the implementation here is highly specialized to
%     Jacobian objecs (e.g., MT2Dsens);  returns a ProjSensMatrix object
%
%   Gary D. Egbert, 2010
%   College of Oceanic and Atmospheric Sciences
%   Oregon State University
%   make sure input Jacobian is normalized
    if nargin ==4
       saveMTX = false;
    end
    
    if ~J.normalized
        Normalize(J,CmHalf);
    end
   if ~dHat.normalized
        dTilde = Normalize(dHat,1);
   end
    
    %   create projected sensitivity matrix object
    Jproj = ProjSensMTX(dTilde.NTX,K,saveMTX);
    Jproj.d = dHat;
    
    beta_alpha = zeros(3,1);
    %  Bi-diag1 of Page and Saunders
    b = dTilde;
    beta = sqrt(dot(b,b));
    beta_alpha(3)  = abs(beta);
    U = (1./beta)*b;
    if saveMTX
        mMTX = JT_times_d_MTX(J,U);
        m = mMTX{1};
        for k = 2:length(mMTX)
            m = m+mMTX{k};
        end
    else
       m = JT_times_d(J,U);
    end
    alpha = sqrt(dot(m,m));
    V = (1./alpha)*m;
    beta_alpha(1) = beta;
    beta_alpha(2) = alpha;
    %   projected sensitivity matrix saves sensitivities for each
    %   transmitter  separately ... sum  these, and put into a model param
    %   object
    if saveMTX
        UpdateProjSensMTX(Jproj,U,mMTX,beta_alpha);
    else
        UpdateProjSensMTX(Jproj,U,m,beta_alpha);
    end
    
    for k = 2:K
        d = J_times_m(J,V);
        b = linComb(1,d,-alpha,U);
        beta = sqrt(dot(b,b));
        beta_alpha(3) = abs(alpha*beta);
        U = (1./beta)*b;    
        if saveMTX
            mMTX = JT_times_d_MTX(J,U);
            m = mMTX{1};
            for l = 2:length(mMTX)
                m = m+mMTX{l};
            end
        else
            m = JT_times_d(J,U);
        end
        V = linComb(1,m,-beta,V);
        alpha = sqrt(dot(V,V));
        V = (1./alpha)*V;
        beta_alpha(1) = beta;  
        beta_alpha(2) = alpha;
        if saveMTX
            UpdateProjSensMTX(Jproj,U,mMTX,beta_alpha);
        else
            UpdateProjSensMTX(Jproj,U,m,beta_alpha);
        end
    end

    end