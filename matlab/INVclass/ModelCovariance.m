classdef ModelCovariance < handle
        % ABSTRACT base class for model covariances
    properties 
        grid;
    end

    methods (Abstract)
        
         %******************************************************************
         mOut = CovMult(Cm,mIn,nTimes)
         %   apply covariance smoother, computing mOut = Cm * mIn
         %   if nTimes = 2, the operator is applied twice
         
    end    % methods
end   % classdef
