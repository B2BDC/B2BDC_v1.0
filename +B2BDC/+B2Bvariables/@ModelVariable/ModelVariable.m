classdef ModelVariable
   % The object represents a single parameter
   
   % Created: June 17, 2015
   
   properties
      NominalValue  % The nominal value of the parameter
      LowerBound    % The lower bound of the parameter
      UpperBound    % The upper bound of the parameter
      Name          % The name of the parameter
   end
   
   methods
      function obj = ModelVariable(name,LB,UB,val)
         %   obj = ModelVariable(name, LB, UB, val) creates a 
         %   B2BDC.B2Bvariables.ModelVariable object with the properties
         %   Name, LowerBound, UpperBound, and NominalValue specified by the 
         %   input variables name, LB, UB, and val, respectively.
         
         %   obj = ModelVariable(name, LB, UB) creates a
         %   B2BDC.B2Bvariables.ModelVariable object with the properties
         %   Name, LowerBound, and UpperBound specified by the 
         %   input variables name, LB, and UB, respectively. The property
         %   NominalValue is calculated by the average of LB and UB.
         if nargin > 0
            if ischar(name)
               obj.Name = name;
            else
               error('Variable name should be a string');
            end
            if ~isscalar(UB) || ~isscalar(LB)
               error('Uncertainty bound should be scalar')
            elseif LB >= UB
               error(['The domain of variable ' name ' is empty']) 
            else
                obj.LowerBound = LB;
                obj.UpperBound = UB;
            end
            if nargin > 3
               obj.NominalValue = val;
            else
               obj.NominalValue = 0.5*(LB+UB);
            end
         end
      end
   end
   
end

