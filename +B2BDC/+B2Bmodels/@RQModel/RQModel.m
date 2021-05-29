classdef RQModel < B2BDC.B2Bmodels.Model
  % A RQModel is a rational quadratic model object inherited from
  % the B2BDC.B2Bmodel.Model superclass
   
   % Created: June 17, 2015     Wenyu Li
   % Modified: Nov 1, 2015     Wenyu Li (error stats added)
   
   properties
      Numerator = [];  % The symmetric coefficient matrix of the quadratic numerator of size (nVar+1)-by-(nVar+1)
      Denominator = []; % The symmetric coefficient matrix of the quadratic denominator of size (nVar+1)-by-(nVar+1)
      K = []; % Conditional parameter defines the upper bound of denominator over variable domain
   end
     
   methods
      function obj = RQModel(Npoly,Dpoly,vars,k,yScale,err)
         % Constructor of the B2BDC.B2Bmodels.RQModel object.
         %
         % The input arguments are:
         %   Npoly - The coefficient matrix of the quadratic numerator
         %   Dpoly - The coefficient matrix of the quadratic denominator
         %   vars - VariableList object containing all model variables
         %   k - The condition parameter defining the upper bound of the denominator over variable ranges
         %   err - Mean and standard deviation of the fitting error of the RQModel
         
         if nargin > 0
            if ~issymmetric(Npoly) || ~issymmetric(Dpoly)
               error('Numerator matrix and denominator matrix must be symmetric')
            else
               obj.Numerator = Npoly;
               obj.Denominator = Dpoly;
            end 
         end
         if nargin > 2
            obj.Variables = vars;
         end
         if nargin > 3
            if isscalar(k) && k > 1
               obj.K = k;
            else
               error('Parameter K must be a scalar greater than 1')
            end
         end
         if nargin > 4
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
         %   Y = EVAL(OBJ, X) evaluates a rational quadratic model at X
         %   sampled points to produce a column vector Y of the
         %   model's output. X is a matrix of size nSample-by-nVariable,
         %   where nSample is the number of samples to be evaluated and
         %   nVariable is the number of variables in the model.
         %
         %   Y = EVAL(OBJ, X, VAROBJ) also evaluates the rational quadratic
         %   model at X sampled points, where X can be of size greater than
         %   nVariable. VAROBJ will specify which columns of X to be used
         %   in evaluating the model to produce a column vector Y.
         if nargin > 2
            oldVar = obj.VarNames;
            newVar = {varObj.Values.Name}';
            [~,~,id] = intersect(oldVar, newVar, 'stable');
            X = X(:,id);
         end
         if size(obj.Numerator) ~= size(X,2)+1
            error('Wrong input dimension of variables')
         else
            x1 = [ones(size(X,1),1), X];
            xNew = B2BDC.Fitting.expandBasis(x1(:,2:end));
            Nvec = B2BDC.Fitting.coef2vec(obj.Numerator);
            Dvec = B2BDC.Fitting.coef2vec(obj.Denominator);
            num = xNew * Nvec;
            den = xNew * Dvec;
            y = num./den;
         end
      end     
   end
   
end

