classdef DataVector
    % Abstract base class for data vectors; this just defines methods and
    % interfaces that need to be implemented.   The multi-Tx structure is
    % also built into all instances, since this aspect is used at higher
    % levels also

   properties
      %   the properties are not abstract, although the structure of the
      %   data vector for a single transmitter (d) might be anything
      NTX  %    Integer giving number of transmitters
      Nd   %    Integer: total number of real data
      d    %    cell array, one cell for each transmitter
      normalized = false   %   logical variable: is data vector normalized
         %     by standard deviations
   end

   methods (Abstract)
       %**********************************************************************
       [obj] = readVec(obj,cfile)
       %
       %  Usage: [obj,status] = readVec(obj,cfile)
       %
       %   reads data vector from impedance file, puts contents into
       %     empty (but already created) data vector object
 
       %*******************************************************************
       [status] = writeVec(obj,cfile,info,units,isign)
       %  Usage:  [status] = writeZ_2D(cfile,allData,info,units,isign);
       %   arguments info, units and isign are optional.
       %   write contents of cell array allData to file
       %   cfile.  There is one cell per period; each
       %   cell contains all information necessary to define
       %   data (locations, values, error standard dev) for
       %   each period (transmitter)

       %**************************************************************************
      [dOut,dErr,normalized] = ExtractVec(obj)
       %  Makes a standard real vector out of an MT2DZ object
       %  Also returns data standard error as a vector
       %
       %  Usage: [dOut,dErr] = ExtractVec(obj);
       %         [dOut,dErr,normalized] = ExtractVec(obj);
       %    optional third argument returns the normalized attribute stored
       %    in the input data vector

       %**********************************************************************
       objOut = SetVec(objIn,dIn,normalized)
       %  Usage: objOut = SetVec(objIn,dIn);
       %  copies the contents of an simple vector object into an existing MT2DZ
       %    object, keeping sizes, order, station coordinates the same (size of
       %    input vector and MT2DZ object must be compatible)
       %   optional argument "normalized" is used to set attribute in
       %   output vector ... does not actually result in any 
       %   normalization/unnormalization inside routine
  
       %*******************************************************************
       objZeroed = zeroVec(obj)  
       %  zeros the data values in an existing data vector object
       %*******************************************************************
       objSum = plus(obj1,obj2)     
       %  Implements addition for data vector objects ... simple at present,
       %    just copies metadata (and error bars) from first object
             %   so far error bars  are just copied from d1;
       %
       %  Usage : d = plus(d1,d2);
       %          d = d1+d2;
       %*******************************************************************
       objPartialSum = add1Tx(obj,j,c,d1)     
       %  adds c*d1, partial data vector for a single transmitter to data
       %     vector object
       %  Usage : obj = add1Tx(obj,k,c,d1) ;
       %*******************************************************************
       objXc1 = mult1Tx(obj,j,c)     
       %  multiplies data vector for a single transmitter (j) by a constant
       %  c
       %  Usage : obj = mult1Tx(obj,k,c) ;
       %*******************************************************************
       objDiff = minus(obj1,obj2)     
       %  Implements subtraction for data vector objects ... simple at present,
       %    just copies metadata (and error bars) from first object
       %   so far error bars  are just copied from d1;
       %
       %  Usage : d = plus(d1,d2);
       %          d = d1+d2;
       %*******************************************************************
       objXc = mtimes(c,obj1)     
       %  Implements scalar multiplication for data vector objects ... simple at present,
       %    c must be double, must be on the left in the expression obj = c*obj1
       %
       %  Usage : d = c*d1
       %          d = mtimes(c,d1);
       %*******************************************************************
       objlinComb = linComb(c1,obj1,c2,obj2)
       %  Usage : [d] = linCombDat(d1,c1,d2,c2);
       %   computes d = c1*d1+c2*d2 for two data vector object
       %    so far error bars  are just copied from d1;
       %*******************************************************************
       n = lengthDat(obj)
       %  Total length of a data vector object
       %
       %  Usage [n] = lengthDat(d);
       %******************************************************************
       objNormalized = Normalize(obj,Ntimes)
       %  Normalizes a data vector by dividing by error standard deviation
       %      Ntimes is 1 to divide by standard deviation
       %                2 to divide by variance
       %             This is optional, default is 1
       %   NO MODIFICATION TO ERROR STANDARD DEVIATION STORED IN OBJECT
       %******************************************************************
       objUnNormalized = UnNormalize(obj,Ntimes)
       %  "UnNormalizes" a data vector by multiplying by error standard deviation
       %      Ntimes is 1 to multiply by standard deviation
       %                2 to multiply by variance
       %             This is optional, default is 1
       %   NO MODIFICATION TO ERROR STANDARD DEVIATION STORED IN OBJECT
       %******************************************************************
       CdInv = InvDataErr(obj)
       % Makes a standard real vector containing the inverse of the data
       % error standard deviation from an input impedance data vector object
       %
       % Usage:  CdInv = InvDataErr(dIn)
       %*******************************************************************
       objNormalized = MultCdInv(obj,objErr)
       %  Normalizes a data vector by dividing by error standard deviation
       %   whether the object is already normalized or not
       %   if a second argument is present errors from this data vector are
       %   used for the normalization, and these errors are copied into
       %    the output vector
       %*******************************************************************
       ip = dot(obj1,obj2)
       %   REAL inner product of two MT2DZ (possibly complex)
       %    data objects, using error standard deviation 
       %    defined from d1, if both objects are unnormalized 
       %  Usage :  ip = dot(d1,d2);
       %*******************************************************************
       ip = dotCmplx(obj1,obj2)
       %   inner product of two (complex impedance)
       %    data objects, using error standard deviation 
       %    defined from d1 
       %  Usage :  ip = dotCmplx(d1,d2)

       %*******************************************************************
       ip = dotMTX(obj1,obj2)
       %   REAL inner product of two (possibly complex impedance)
       %    data objects, using error standard deviation 
       %    defined from d1 (if both objects are unnormalized  ... 
       %    operates transmitter by transmitter, returns an array in ip
       %  Usage :  ip = dotMTX(d1,d2);
       %*******************************************************************
       ip = dotCmplxNoCov(obj1,obj2)
       %   inner product of two (complex impedance)
       %    data objects, No  normalization by standard deviations
       %  Usage :  ip = ipDat(d1,d2);
       %**************************************************************************
       ip = dotNoCov(obj1,obj2)
       %   REAL inner product of two (possibly complex impedance)
       %    data objects, No normalization by standard deviations
       %  Usage :  ip = ipDat(d1,d2);
       %**************************************************************************
       objUnit = UnitVec(objIn,MTX)
       %   normalizes data vector so that it has unit norm; if optional
       %   argument MTX is present and MTX = true each transmitter is
       %   normalized separately (to have unit norm).
    end    %    methods
end     %   classdef
