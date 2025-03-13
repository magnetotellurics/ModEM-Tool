
classdef Sensitivity < handle
    % Abstract class for sensitivity matrix (Jacobian)
    
    properties
       
        N      %    number of rows ... dimension of real data space
        M      %    number of columns  .... dimension of model space
        d      %   data vector template; length (considered as real vector)
               %     will be same as total number of sensitivities, and
               %     sensitivities are stored in same order;  allows us to
               %     sort out what is what in the simple cell array of
               %     sensitivities.   Better structures in the future!
               %     ALSO: keep original data errors here!
        grid   %   This is just a handle of the grid object, making it
               %     easier to plot sensitivities, etc.
        paramType   %   paramter type (all sensitivities would be the same)
        normalized = false
        dataVecClass  %    data vector class for this
                                 %    implementation of sensitivity
        modelVecClass   %  model parameter class
        J_stored   %   is full sensitivity stored or not ... so far just
                  %   for information
    end    
        
    methods (Abstract)
       
        %******************************************************************
        updateSens(obj,m)
        %    upate sensitivity object for model parameter m
        %
        %    Usage: updateSens(obj,m)
        %
        %    Assumes object has been created and data vector structure has
        %          been set
        %******************************************************************
        dTilde = ExtractNormVec(obj,dIn)
        %******************************************************************
        b = SetNormVec(obj,bTilde)
        %******************************************************************
        [dOut] = J_times_m(obj,m)
        %  Usage : [dOut] = J_times_m(J,m);
        %    multiply m by sensitivity matrix, put result in data space
        %    object dOut

        %******************************************************************
        [mOut] = JT_times_d(obj,d)
        %  Usage : [mOut] = JT_times_d(J,d);
        % multiply by transpose of sensitivity matrix (stored as
        %     cell array of data/model space structures)

        %******************************************************************
        m =  rowJ(obj,k)
        %   extract row k from representer
 
        %******************************************************************
        obj2 = copy(obj1)
        %   makes a copy of actual sensitivity object (not handle)

        %******************************************************************
        Normalize(obj,CmHalf)
        %   Normalizes Jacobian structure, overwritting input Jacobian by
        %   J_tilde = C_m^{1/2} J C_d^{-1/2} and d_tilde = C_d{-1/2} d

        %******************************************************************    
        [g] = JT_times_d_MTX(obj,d,SaveCmplx)
        %  variant on JT_times_d ... does not sum over
        %   frequencies, and optionally returns "imaginary part",
        %   for transmitters in input data vector d which are complex
        %   output is a cell array of size NTX or (NTX,2) for real and 
        %   complex cases, respectively
 
    end    %  methods ....
        
end    %   classdef
