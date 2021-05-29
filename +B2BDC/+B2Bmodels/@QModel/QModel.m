classdef QModel < B2BDC.B2Bmodels.Model
   % A QModel is a quadratic model object inherited from
   % B2BDC.B2Bmodel.Model superclass
   
   % Created: June 17, 2015     Wenyu Li
   % Modified: Nov 1, 2015     Wenyu Li (error stats added)
   properties 
      CoefMatrix = [];  % The original symmetric coefficient matrix of the quadratic model of size (nVar+1)-by-(nVar+1)
   end
   
   properties (Dependent, Hidden = true)
      Hessian;
   end
   
   methods
      function obj = QModel(coef,vars,err)
            % Constructor of the B2BDC.B2Bmodels.QModel object.
            %
            % The input arguments are:
            %   coef - The coefficient matrix of the quadratic model
            %   vars - VariableList object containing all model variables
            %   err - Mean and standard deviation of the fitting error of the QModel
         if nargin > 0
            if issymmetric(coef)
               obj.CoefMatrix = coef;
            else
               obj.CoefMatrix = 0.5*(coef+coef');
               warning('Coefficient matrix must be symmetric')
            end
         end
         if nargin > 1
            if size(obj.CoefMatrix,1) == vars.Length+1
               obj.Variables = vars;
            else
               error('Mismatched coefficient matrix and variable dimension');
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
        %   Y = EVAL(OBJ, X) evaluates a quadratic model at X sampled
        %   points to produce a column vector Y of the model's output. X
        %   is a matrix of size nSample-by-nVariable, where nSample is the
        %   number of samples to be evaluated and nVariable is the number
        %   of variables in the model.
        %
        %   Y = EVAL(OBJ, X, VAROBJ) also evaluates the quadratic model at
        %   X sampled points, where X can be of size greater than nVariable.
        %   VAROBJ will specify which columns of X to be used in evaluating
        %   the model to produce a column vector Y of the model's output.
        
         if nargin > 2
            oldVar = obj.VarNames;
            newVar = {varObj.Values.Name}';
            [~,~,id] = intersect(oldVar, newVar, 'stable');
            X = X(:,id);
         end
         if size(obj.CoefMatrix,1) ~= size(X,2)+1
            error('Wrong input dimension of variables')
         else
%             nSample = size(X,1);
%             x1 = [ones(nSample,1), X];
            xNew = B2BDC.Fitting.expandBasis(X);
            covec = B2BDC.Fitting.coef2vec(obj.CoefMatrix);
            y = xNew*covec;
%             coef = obj.CoefMatrix;
%             y = diag(x1 * coef * x1');
         end
      end
      
      function y = get.Hessian(obj)
         y = 2 * obj.CoefMatrix(2:end,2:end);
      end
   end
   
   methods (Hidden = true)
      function [A,b] = calGradient(obj)
         quadCoef = obj.CoefMatrix(2:end,2:end);
         A = 2 * quadCoef;
         b = 2 * obj.CoefMatrix(2:end,1);
      end
   end
   
end

