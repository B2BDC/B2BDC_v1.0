classdef SampleOption < handle
   % A B2BDC.Option.SampleOption object that specifies the sampling crieteria.
   % The syntax of the input and the meaning of each property are in the
   % following:
   % opt = B2BDC.Option.SampleOption(PropName, PropValue ... )
   % ---------------------------------------------------------------------
   % SampleMethod:
   %   'Gibbs' - Gibbs sampler
   %   'HR' - Hit-and-Run sampler
   %   'AM' - Adaptive Metropolis sampler
   %   'CWAM' - Component-wise Adaptive Metropolis sampler
   % ---------------------------------------------------------------------
   % StepInterval:
   %   a numerical scalar specifies the thinning parameter, i.e., the number of
   %   samples that are discarded between two returning sampled points plus one.
   % ---------------------------------------------------------------------
   % AMParam:
   %   a structure variable contains information associated with the
   %   Adaptive Metropolis and Component-wise Adaptive Metropolis samplers
   %      WarmStart - a structure array specifies the recorded sampling information related to
   %                  the AM and CWAM sampler
   %      NumOfSigma - a scalar indicates how many standard deviations
   %                   correspond to a half uncertainty interval length
   %      Sigma -  a nVar-by-nVar matrix specifies the starting covariance
   %               matrix for generating the candidate directions
   % ---------------------------------------------------------------------
   % Created: Oct 11, 2016    Wenyu Li
   %  Modified: Septempber 2, 2015   Wenyu Li (Data normalization option added)
   
   
   properties
      SampleMethod = [];
      StepInterval = [];
      AMparam = [];
   end
   
   methods
      function  obj = SampleOption(inputCell)
         % To generate a B2BDC.Option object
         p = {'SampleMethod','StepInterval','AMparam'};
         if nargin > 0
            nin = length(inputCell);
         else
            nin = 0;
         end
         if mod(nin,2) ~= 0
            error('Wrong number of input argument')
         end
         nset = floor(0.5*nin);
         if nset > 0
            for i = 1:nset
               if any(strcmp(p,inputCell{2*i-1}))
                  id = find(strcmp(p,inputCell{2*i-1}) == 1);
                  switch id
                     case 1
                        obj.SampleMethod = inputCell{2*i};
                     case 2
                        obj.StepInterval = inputCell{2*i};
                     case 3
                        obj.AMparam = inputCell{2*i};
                  end
               else
                  error('Invalid input property names')
               end
            end
         end
         if isempty(obj.SampleMethod)
            obj.SampleMethod = 'Gibbs';
         end
         if isempty(obj.StepInterval)
            obj.StepInterval = 1;
         end
         if isempty(obj.AMparam)
            cc.WarmStart = [];
            cc.NumOfSigma = 2;
            cc.Sigma = [];
            obj.AMparam = cc;
         end
      end
      
      function set.SampleMethod(obj,str1)
         c = {'Gibbs','HR','AM','CWAM'};
         if any(strcmp(c,str1))
            obj.SampleMethod = str1;
         else
            error('Invalid input property value')
         end
      end
      
      function set.StepInterval(obj,num1)
         c = round(num1);
         if c > 0
            obj.StepInterval = c;
         else
            error('Invalid input property value')
         end
      end
      
      function set.AMparam(obj,c)
         if isstruct(c) && all(isfield(c,{'WarmStart','NumOfSigma','Sigma'}))
            obj.AMparam = c;
         else
            error('The specified AMparam peroperty has incorrect structure field(s)')
         end
      end
      
   end
   
   methods (Hidden = true)
      % Matlab additional Functions
      y = addlistener() %matlab add listener events
      y = delete() % matlab delete handle
      y = findobj() % matlab find objects with a specified property value
      y = findprop() % matlab find matlab property of handle obj
      y = notify()
      y = le() % matlab less than or equal to
      y = lt() % matlab less than
      y = ne() % matlab not equal to
      y = ge() % matlab greater than or equal to
      y = gt() % matlab greater than
      y = eq() % matlab equal to
      %        y = isvalid() % matlab method for timer obj -- breaks when
      %hidden...
   end
   
end