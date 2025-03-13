classdef ModelParameter
        % ABSTRACT base class fo model parameter vectors; 
        %         this just defines methods and interfaces that need to be 
        %         implemented  in specific instances


    properties 
        %    these properties are inherited by all subclasses
        %   make covariance, grid handles public, so these can accessed
        cov         %   handle of model covariance object   ... hmmm
        grid        %   handle of numerical grid object
    end

    methods (Abstract)
   
       %*******************************************************************
       M = ModelParamLength(obj)
           % returns length of real vector used to represent model
           % parameter
       %*******************************************************************
       M = ModelParamMax(obj)
           %  estimate maximum model parameter size 
       %*******************************************************************
       [obj] = readVec(obj,cfile)
       %   reads model parameter from file cfile
       %*******************************************************************
       [status] = writeVec(obj,cfile)
       %  writes ModelParameter object to file cfile
       %*******************************************************************
       mVec = ExtractVec(obj)
       %   extracts model parameter values and returns as a standard vector
       %     of length ModelParamLength(obj)
        %*******************************************************************
       obj = SetVec(obj,v)
       %   inserts  standard vector v into model parameter obj (reverses 
       %         ExtractVec) 
       %*******************************************************************
       obj = plus(obj1,obj2)
       %    add two model parameters obj = obj1+obj2
       %*******************************************************************
       obj = minus(obj1,obj2)
       %    subtract two model parameters obj = obj1-obj2
       %*******************************************************************
       [m] = linComb(c1,m1,c2,m2)
       %   computes linear combination of two model parameters
       %                 m = c1*m1+c2*m2
       %*******************************************************************
       [mOut] = mtimes(c,mIn)
       % scalar multiplication of a model space object by a real scalar;
       %*******************************************************************
       [ip] = dot(m1,m2)
       %  dot product for (real) model parameters m1, m2
       %*******************************************************************
       [obj] = ZeroVec(obj)
       %    zero existing model parameter 
       %*******************************************************************
       [d] = fwd(m,d0)
       %  Computes predicted data for model parameter m, using d0 as
       %  template
       %***********************************************
       diff  = relDiff(obj1,obj2)
       %   compute relative difference of two model parameter objects,
       %   (obj1-obj2)/sqrt(dot(obj1,obj2))
       
    end    % methods
end   %   classdef
