classdef LModel < B2BDC.B2Bmodels.Model
   % A QModel is a quadratic model object inherited from
   % B2BDC.B2Bmodel.Model superclass
   
   % Created: Sep 20, 2016     Wenyu Li
   
   properties
      CoefVec = []; % The original coefficient vector of the linear model
   end
   
   methods
       function obj = LModel(coef,vars,err)
           % Constructor of the B2BDC.B2Bmodels.LModel object.
           % The input arguments are:
           %   coef - The coefficient vector of the linear model
           %   vars - VariableList object containing all model variables
           %   err - Mean and standard deviation of the fitting error of the QModel
           if nargin > 0
               if size(coef,1) == 1
                   obj.CoefVec = coef';
               else
                   obj.CoefVec = coef;
               end
           end
           if nargin > 1
               if length(obj.CoefVec) ~= vars.Length+1
                   error('Mismatched variable and coefficient dimension');
               else
                   obj.Variables = vars;
               end
           end
           if nargin > 2
               obj.ErrorStats.absMax = err.absMax;
               obj.ErrorStats.absAvg = err.absAvg;
               obj.ErrorStats.relMax = err.relMax;
               obj.ErrorStats.relAvg = err.relAvg;
           else
              obj.ErrorStats.absMax = 0;
              obj.ErrorStats.absAvg = 0;
              obj.ErrorStats.relMax = 0;
              obj.ErrorStats.relAvg = 0;
           end
      end
      
      function y = eval(obj, X, varObj)
         %   Y = EVAL(OBJ, X) evaluates a linear model at X sampled
         %   points to produce a column vector Y of the model's output. X
         %   is a matrix of size nSample-by-nVariable, where nSample is the
         %   number of samples to be evaluated and nVariable is the number
         %   of variables in the model.
         %
         %   Y = EVAL(OBJ, X, VAROBJ) also evaluates the linear model at
         %   X sampled points, where X can be of size greater than nVariable.
         %   VAROBJ will specify which columns of X to be used in evaluating
         %   the model to produce a column vector Y of the model's output.
         if nargin > 2
            oldVar = obj.VarNames;
            newVar = {varObj.Values.Name}';
            [~,~,id] = intersect(oldVar, newVar, 'stable');
            X = X(:,id);
         end
         if length(obj.CoefVec) ~= size(X,2)+1
            error('Wrong input dimension of variables')
         else
            nSample = size(X,1);
            x1 = [ones(nSample,1), X];
            y = x1*obj.CoefVec;
         end
      end
   end
   
end

