classdef Mod3DMTCov < handle
        % another model covariance class for 3DMT ... this one just uses
        % the default ModEM F95 code, calling Mod3DMT to apply the
        % covariance operator  (allows exact comparison with Fortran
        % inversion code)

    properties 
        
    end
    
    methods
        function obj = Mod3DMTCov()
        %   class constructor ... attach grid handle, set default
        %        covariance parameter values, and setup horizontal and
        %        verical smoothing matrices
        end
         %******************************************************************
        function mOut = CovMult(obj,mIn,nTimes)
        %  Usage  : mout = CovMult(sIn,smthParams)
        %  Input : mIn = input conductivity structure
        %          obj = model covariance object
        %          nTimes = optional argument: number of times to apply 
        %		           the smoother (default is 1)
        %  Output : sOut = smoothed output conductivity structure
            if nargin == 2
                nTimes = 1;
            end
            %    make scratch directory if it doesn't exist ...
             if exist('scratch','dir')==0
                 mkdir('scratch');
             end
             %   write scratch file for input conductivity parameter
             m0_File = 'scratch/m0.cpr';
             writeVec(mIn,m0_File);
             %   file for smoothed model parameter
             m1_File = 'scratch/m0.cpr';
             %   call Mod3DMT to apply covariance smoother
             status = system(['Mod3DMT -C FWD ' m0_File ' ' m1_File]);
             [mOut] = readVec(mIn,m1_File);
             if nTimes == 2
                 %  repeat
                 writeVec(mOut,m0_File);
                 status = system(['Mod3DMT -C FWD ' m0_File ' ' m1_File]);
                 [mOut] = readVec(mIn,m1_File);
             end    
        end
    end    % methodss
end   % classdef
