classdef InvIterControl

   properties
      %   properties used by a wide range of iterative schemes
      %   below are defaults, used if no arguments are provided to constructor
      niter = 0;
      iterMax = 20;
      RMStolerance = 1.05;
      ResidualTolerance = .01; 
      relErr ;         %   array (indexed by iteration #) of relative error
      RMS              %   array (indexed by iteration #) of RMS misfit

      %   additional properties for specialized uses .... 

          %    Occam style inversions
          NU               %   array of tradeoff parameters
          X2               %   array of squared misfits
          OccamPhase       %   phase I or phase II
          PhaseI_Max = 10   %   max # of Occam phase I iterations
          PhaseII_Max = 1  %   max # of Occam phase II iterations
          ModifiedOccam=false   %   true to do modified Occam, relaxing
                           %   elimination of model  component not in span
                           %   of current representers
         X2mu
         MU

          %    NLCG inversions
          PEN              %   array (indexed by iteration #) of Penalty
                           %         functional value
          penChange    
          nCGmax = 10
          minChange = .00001
          UserData
   end

   methods 
   % class constructor
   function obj = InvIterControl(varargin)
      %  constructor for class IterControl
      if nargin > 0
         n = length(varargin);
         if mod(n,2)
            error(1,'%s\n','Optional arguments to IterControl must occur in pairs')
         end
         for k = 1:2:n
            option = lower(varargin{k});
            switch option
               case 'itermax'
                  obj.iterMax = varargin{k+1};
               case 'tolerance'
                  obj.ResidualTolerance = varargin{k+1};
               otherwise
                  error(1,'%s\n','Optional argument to IterControl not defined')
            end
         end
         obj.niter = 0;
         obj.relErr = zeros(obj.iterMax,1);
         obj.UserData = cell(obj.iterMax,1);
      end
      %   else just use defaults
   end
   %***********************************************************************
   function DONE = checkConvergence(obj)

      DONE = ( obj.relErr(obj.niter) < obj.ResidualTolerance ) || ...
                 (obj.niter > obj.iterMax);
   end
   %***********************************************************************
   function DONE = FitToTolerance(obj)

      DONE = ( obj.RMS(obj.niter) < obj.RMStolerance ) || ...
                    (obj.niter > obj.iterMax);
   end
   %***********************************************************************
   function DONE = NLCGconverged(obj)

      DONE = ( obj.RMS(obj.niter) < obj.RMStolerance ) || ...
                    (obj.niter > obj.iterMax) || ...
                    (obj.penChange < obj.minChange);
   end
   %***********************************************************************
   function nTot = niterTotal(obj)
       nOuter = length(obj);
       nTot = 0;
       for k = 1:nOuter
           nTot = nTot+obj(k).niter-1;
       end
   end
   %***********************************************************************
   function hfig = plotOccCurves(obj)
       nOuter = length(obj);
       hfig = figure('PaperPosition',[1,1,5,5],'Position',[100,100,500,500]);
       %lns = {'-','--','-','--','-','--','-','--'};
       %clrs = [0,0,.5;0,0,.5;.5,0,0;.5,0,0;0,0,0;0,0,0;0.,.5,0;0,.5,0];
       lns = {'-','-','-','-','-','-','-','-'};
       clrs = [0; .2 ; .4 ; .5 ; .6 ; .7 ; .75 ; .8]*ones(1,3);
  
       for k = 1:nOuter
          k1 = min(length(lns),k);
          NUsm = smooth(obj(k).NU,1);
          X2sm = smooth(obj(k).X2,1);
          loglog(NUsm,sqrt(X2sm)/1.4,'color',clrs(k1,:),...
              'linewidth',3,'linestyle',lns{k1});
          hold on
       end
       set(gca,'fontweight','demi','fontsize',14)
   end
   %***********************************************************************
   function hfig = plotCvgc(obj)
       nOuter = length(obj);
       hfig = figure('PaperPosition',[1,1,5,5],'Position',[100,100,500,500]);
       lns = {'-','--','-','--','-','--','-','--'};
       clrs = [0,0,.7;0,0,.7;.7,0,0;.7,0,0;0,0,0;0,0,0;0.,.7,0;0,.7,0];
       lns = {'-','-','-','-','-','-','-','-'};
       
       clrs = [0; .2 ; .4 ; .5 ; .6 ; .7 ; .75 ; .8]*ones(1,3);
       for k = 1:nOuter
          k1 = min(length(lns),k);
          semilogy([1:obj(k).niter],obj(k).relErr,'color',clrs(k1,:),...
              'linewidth',3,'linestyle',lns{k1});
          hold on
       end
       set(gca,'fontweight','demi','fontsize',14)
   end      
       
       
   end   % methods

end   % class

