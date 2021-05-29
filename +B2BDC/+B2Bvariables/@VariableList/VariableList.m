classdef VariableList < B2BDC.Util.IContainer
% This object represents a list of parameters
   
% Created: June 17, 2015, myf + Wenyu Li
% Modified: Nov 1, 2015     Wenyu Li (Sample function added)
   
   properties (SetAccess = private)
      ExtraLinConstraint = struct('A',[],'LB',[],'UB',[]);  % save the information of linear constraints of the parameter vector in addition to the specified interval uncertainty of each individual parameter
   end
   
   methods
      function obj = VariableList()
         %   obj = VariableList() creates a B2BDC.B2Bvariables.VariableList
         %   object, which acts as a container for the ModelVariable 
         %   objects.
         valClass = 'B2BDC.B2Bvariables.ModelVariable';
         obj = obj@B2BDC.Util.IContainer(valClass);
      end
      
      function obj = addList(obj,varList)
         %   newobj = addList(obj,varList) adds the variable list saved in varList to the variable list
         %   saved in obj and returns the extended variable list as newobj. If any parameters are
         %   in common between obj and varList, the intersection of their respective uncertainties
         %   is used in the extended parameter list.
         A1 = obj.ExtraLinConstraint.A;
         LB1 = obj.ExtraLinConstraint.LB;
         UB1 = obj.ExtraLinConstraint.UB;
         if isempty(obj.Values)
            name1 = {};
         else
            name1 = {obj.Values.Name};
         end
         A2 = varList.ExtraLinConstraint.A;
         LB2 = varList.ExtraLinConstraint.LB;
         UB2 = varList.ExtraLinConstraint.UB;
         name2 = {varList.Values.Name};
         obj = obj.clearExtraConstraint;
         variables = varList.Values;
         for i = 1:varList.Length
            variable = variables(i);
            [y,ind] = obj.find('Name',variable.Name);
            if isempty(y)
               obj = obj.add(variable);
            else
               if variable.LowerBound > obj.Values(ind).LowerBound
                  obj.Container{ind}.LowerBound = variable.LowerBound;
               end
               if variable.UpperBound < obj.Values(ind).UpperBound
                  obj.Container{ind}.UpperBound = variable.UpperBound;
               end
            end
         end
         allName = {obj.Values.Name};
         nV = obj.Length;
         if ~isempty(A1)
            n1 = size(A1,1);
            At = zeros(n1,nV);
            [~,~,id] = intersect(name1,allName,'stable');
            At(:,id) = A1;
            obj = obj.addLinearConstraint(At,LB1,UB1);
%             obj = obj.addLinearConstraint([At;-At],[UB1;-LB1]);
         end
         if ~isempty(A2)
            n2 = size(A2,1);
            At = zeros(n2,nV);
            [~,~,id] = intersect(name2,allName,'stable');
            At(:,id) = A2;
            obj = obj.addLinearConstraint(At,LB2,UB2);
%             obj = obj.addLinearConstraint([At;-At],[UB2;-LB2]);
         end
         if ~obj.checkFeasibility
            disp('The resulted variableList is not feasible')
         end
      end
      
      function xSample = makeDesignSample(obj,nSample,method)
         %   xSample = makeDesignSample(obj, nSample, method) generates nSample i.i.d (independently and identically
         %   distributed) samples using the approach specified by method, ('lhs', 'sobol', or 'uniform', case insensitive).
         
         %   xSample = makeDesignSample(obj, nSample) generates nSample i.i.d (independently and identically
         %   distributed) samples using the 'lhs' (Latin-hypercube sampling) approach.
         if nargin < 3
             xSample = makeLHSsample(obj,nSample);
         else
             switch(lower(method))
                 case 'lhs'
                     xSample = makeLHSsample(obj,nSample);
                 case 'sobol'
                     xSample = makeSobolsample(obj,nSample);
                 case 'uniform'
                     xSample = makeUniformsample(obj,nSample);
                 otherwise
                     error('Misspecified sampling method');
             end
         end
      end
      
      function obj = deleteVariable(obj,varIdx)
         %   newobj = deleteVariable(obj, varIdx) removes the variable specified by varIdx in the
         %   variable list specified by obj. The new variable list is returned as newobj.
         if ischar(varIdx)
            [~,id,~] = intersect(varIdx,{obj.Values.Name});
         elseif isscalar(varIdx)
            id = varIdx;
         else
            error('Wrong input index type or length')
         end
         A = obj.ExtraLinConstraint.A;
         tA = A(:,id);
         if any(tA)
            error('Extra linear constraints contain the variables you want to delete')
         end
         if ~isempty(id)
            obj = obj.remove(id);
         else
            disp('The variable you want to delete is not in the current list')
         end
      end
      
      function varBD = calBound(obj)
         %   varBD = calBound(obj) returns the prior interval uncertainty of each individual variable as a 2-D array varBD. 
         %   The first and second columns of varBD are the lower and upper bounds, respectively.
         varVal = obj.Values;
         varBD = [[varVal.LowerBound]', [varVal.UpperBound]'];
      end
      
      function varOB = calNominal(obj)
         %   varOB = calBound(obj) returns the nominal value of each individual variable as a 1-D
         %   column array varOB.
         varOB = [obj.Values.NominalValue]';
      end
   
      function obj = changeBound(obj, newBD, idx)
         %   newobj = changeBound(obj, newBD, idx) updates the uncertainty bounds (if newBD
         %   has two columns) and the nominal values (if newBD has three columns and the last
         %   column corresponds to the new nominal value) in obj and returns the updated variable
         %   list as newobj. The variables updated are specified by the input idx, a cell array of the
         %   variable names or a numerical array of the positions of the variables in the list.
         
         %   newobj = changeBound(obj, newBD) updates the uncertainty bounds (if newBD has
         %   two columns) and the nominal values (if newBD has three columns and the last column
         %   corresponds to the new nominal value) of all variables in obj and returns the updated
         %   variable list as newobj.
         allVarName = {obj.Values.Name};
         if nargin == 3
            if length(idx) ~= size(newBD,1)
               error('Inconsistent new bounds size with number of variables')
            elseif iscell(idx)
               for i = 1:length(idx)
                  [~,id] = intersect(allVarName, idx{i});
                  if isempty(id)
                     error(['The variable ' idx{i} ' is not in the VariableList'])
                  end
                  varName = obj.Values(id).Name;
                  if size(newBD,2) == 2
                     newVar = B2BDC.B2Bvariables.ModelVariable(varName, newBD(i,1), newBD(i,2));
                  else
                     newVar = B2BDC.B2Bvariables.ModelVariable(varName, newBD(i,1), newBD(i,2), newBD(i,3));
                  end
                  obj = obj.replace('Name',varName,newVar);
               end
            elseif isnumeric(idx)
               for i = 1:length(idx)
                  varName = obj.Values(idx(i)).Name;
                  if size(newBD,2) == 2
                     newVar = B2BDC.B2Bvariables.ModelVariable(varName, newBD(i,1), newBD(i,2));
                  else
                     newVar = B2BDC.B2Bvariables.ModelVariable(varName, newBD(i,1), newBD(i,2), newBD(i,3));
                  end
                  obj = obj.replace('Name',varName,newVar);
               end
            else
               error('Wrong input index')
            end
         elseif nargin == 2 && size(newBD,1) == obj.Length
            varName = {obj.Values.Name};
            if size(newBD,2) == 2
               Var = generateVar(varName,newBD);
            else
               Var = generateVar(varName,newBD(:,1:2),newBD(:,3));
            end
            if ~isempty(obj.ExtraLinConstraint.A)
               A = obj.ExtraLinConstraint.A;
               LB = obj.ExtraLinConstraint.LB;
               UB = obj.ExtraLinConstraint.UB;
               obj = Var.addLinearConstraint(A,LB,UB);
            end
         else
            error('Wrong number of input arguments')
         end
      end
      
      function newVar = addLinearConstraint(obj,A,b,c)
         %   newVar = addLinearConstraint(obj,A,b) adds linear constraints to the variables such
         %   that Ax <= b. Since the constraints are saved as LB <= Ax <= UB in the object, the
         %   corresponding lower bounds are calculated by solving a linear programming.
         
         %    newVar = addLinearConstraint(obj,A,b,c) adds linear constraints to the variables such
         %    that b <= Ax <= c.
         if nargin < 4
            if size(b,1) == 1
               b = b';
            end
            newVar = obj.clearExtraConstraint;
            [nC, nVar] = size(A);
            if nVar ~= obj.Length
               error('Invalid input matrix dimension')
            end
            if nC ~= length(b)
               error('Invalid input vector dimension')
            end
            [idA,id] = intersect(A,-A,'rows','stable');
            nint = 0.5*size(idA,1);
            nuni = nC - 2*nint;
            A0 = zeros(nint+nuni,nVar);
            LB = zeros(nint+nuni,1);
            UB = zeros(nint+nuni,1);
            if ~isempty(idA)s
               count = 1;
               for i = 1:2*nint
                  idx = find(idA(i,:),1);
                  if idA(i,idx) > 0
                     A0(count,:) = idA(i,:);
                     [~,id1] = intersect(A,idA(i,:),'rows');
                     [~,id2] = intersect(-A,idA(i,:),'rows');
                     UB(count) = b(id1);
                     LB(count) = -b(id2);
                     count = count+1;
                  end
               end
            end
            tmpA = A;
            tmpb = b;
            A(id,:) = [];
            b(id,:) = [];
            A0(nint+1:end,:) = A;
            UB(nint+1:end) = b;
            for i = 1:nuni
               tmpC = A(i,:);
               LB(nint+i) = finddirectionmin(obj,tmpA,tmpb,tmpC);
            end
            for i = 1:nuni
               idx = find(A0(nint+i,:),1);
               if A0(nint+i,idx) < 0
                  A0(nint+i,:) = -A(i,:);
                  tmpLB = -UB(nint+i);
                  UB(nint+i) = -LB(nint+i);
                  LB(nint+i) = tmpLB;
               end
            end
            if ~isempty(obj.ExtraLinConstraint.A)
               Ai = obj.ExtraLinConstraint.A;
               Li = obj.ExtraLinConstraint.LB;
               Ui = obj.ExtraLinConstraint.UB;
               [~,id1,id2] = intersect(Ai,A0,'rows','stable');
               if ~isempty(id1)
                  %                LB(LB(id2) < Li(id1)) = Li(id1);
                  LB(id2) = max([LB(id2), Li(id1)],[],2);
                  %                UB(UB(id2) > Ui(id1)) = Ui(id1);
                  UB(id2) = min([UB(id2), Ui(id1)],[],2);
               end
               Ai(id1,:) = [];
               Li(id1) = [];
               Ui(id1) = [];
               A0 = [A0;Ai];
               LB = [LB;Li];
               UB = [UB;Ui];
            end
            newVar.ExtraLinConstraint.A = A0;
            newVar.ExtraLinConstraint.LB = LB;
            newVar.ExtraLinConstraint.UB = UB;
         else
            if size(A,1) ~= length(b)
               error('Wrong lower bound dimension')
            end
            if size(A,1) ~= length(c)
               error('Wrong upper bound dimension')
            end
            newVar = obj;
            A0 = obj.ExtraLinConstraint.A;
            b0 = obj.ExtraLinConstraint.LB;
            c0 = obj.ExtraLinConstraint.UB;
            if ~isempty(A0)
               [~,id] = setdiff(A,A0,'rows');
               if ~isempty(id)
                  newVar.ExtraLinConstraint.A = [A0;A(id,:)];
                  newVar.ExtraLinConstraint.LB = [b0;b(id)];
                  newVar.ExtraLinConstraint.UB = [c0;c(id)];
               end
            else
               newVar.ExtraLinConstraint.A = A;
               newVar.ExtraLinConstraint.LB = b;
               newVar.ExtraLinConstraint.UB = c;
            end
         end
         if ~newVar.checkFeasibility
            disp('The resulted variableList is not feasible')
            newVar = obj;
         end
      end
      
      function newVar = makeSubset(obj,varIdx)
         %   newVar = makeSubset(obj,varIdx) generates a new variable list that contains a subset
         %   of the original variable list. The subset is specified by the names contained in the cell
         %   array varIdx or by the index contained in the numerical array varIdx.
         vName = {obj.Values.Name}';
         if ischar(varIdx)
            [~,~,id] = intersect(varIdx,{obj.Values.Name},'stable');
         elseif isvector(varIdx)
            id = varIdx;
         else
            error('Wrong input index type or length')
         end
         H = obj.calBound;
         ob = [obj.Values.NominalValue]';
         newName = vName(id);
         newH = H(id,:);
         newOB = ob(id,:);
         newVar = generateVar(newName, newH, newOB);
         if ~isempty(obj.ExtraLinConstraint.A)
            A = obj.ExtraLinConstraint.A;
            LB = obj.ExtraLinConstraint.LB;
            UB = obj.ExtraLinConstraint.UB;
            A_new = A(:,id);
            idA = true(size(A,1),1);
            for i = 1:size(A,1)
               a1 = sum(A(i,:)~=0);
               a2 = sum(A_new(i,:)~=0);
               if a1 == a2
                  idA(i) = false;
               end
            end
            A_new(idA,:) = [];
            if ~isempty(A)
               LB(idA) = [];
               UB(idA) = [];
               newVar = newVar.addLinearConstraint(A_new,LB,UB);
            end
         end
      end
      
      function newVar = clearExtraConstraint(obj)
         % clear extra constraints
         newVar = obj;
         newVar.ExtraLinConstraint = struct('A',[],'LB',[],'UB',[]);
      end
      
      function xNew = changeCoordinate(obj,xOld,newVList)
         %   xNew = changeCoordinate(obj,xOld,newVList) creates a new 2D numerical array xNew, whose corresponding parameter 
         %   values are assigned from the input xOld for those overlapped between the two variable lists. Parameters not overlapped 
         %   are assigned to zeros.
         name1 = {newVList.Values.Name}';
         name2 = {obj.Values.Name}';
         n1 = newVList.Length;
         n2 = obj.Length;
         if n2 ~= size(xOld,2)
             error('Wrong input sample size')
         end
         ns = size(xOld,1);
         [~,id1,id2] = intersect(name1,name2);
         xNew = zeros(ns,n1);
         if ~isempty(id1)
             xNew(:,id1) = xOld(:,id2);
         end
      end
   end
   
   methods (Static, Hidden = true)
       xval = q2sample(Mgd,idx,H,Xvals,V);
       xVal = quadSample(Q,H,x0,V)
   end
   
   methods (Hidden = true)
      function obj = addVariable(obj,varName,varRange,varValue)
         if nargin > 3
            variable = B2BDC.B2Bvariables.ModelVariable(varName,varRange(1),varRange(2),varValue);
         else
            variable = B2BDC.B2Bvariables.ModelVariable(varName,varRange(1),varRange(2));
         end
         [y,ind] = obj.find('Name',variable.Name);
         if isempty(y)
            obj = obj.add(variable);
         else
            if variable.LowerBound > obj.Values(ind).LowerBound
               obj.Container{ind}.LowerBound = variable.LowerBound;
            end
            if variable.UpperBound < obj.Values(ind).UpperBound
               obj.Container{ind}.UpperBound = variable.UpperBound;
            end
         end
         if ~isempty(obj.ExtraLinConstraint.A)
            obj.ExtraLinConstraint.A = [A, zeros(size(A,1),1)];
         end
      end
      
      function xSample = makeLHSsample(obj,nSample)
         %   XSAMPLE = MAKELHSSAMPLE(OBJ, NSAMPLE) returns latin hypercube 
         %   samples of the VariableList OBJ, within the Variable 
         %   LowerBound and UpperBound. NSAMPLE specifies the number of 
         %   sample points to return from the VariableList domain.
         nVar = obj.Length;
         H = [[obj.Values.LowerBound]', [obj.Values.UpperBound]'];
         dH = diff(H');
         xdesign = lhsdesign(nSample,nVar,'criterion','none');
         xSample = repmat(H(:,1)',nSample,1) + repmat(dH,nSample,1).*xdesign;
      end
      
      function xCand = makeSobolsample(obj,n)
          % XCAND = SOBOLSAMPLE(OBJ,N) returns a n-by-nVar sample matrix. The method
          % generates the sample with sobol sequence.
          nVar = obj.Length;
          sob = sobolset(nVar,'Skip',1e3*randi(1e3),'Leap',1e2*randi(5e2));
          sob = scramble(sob,'MatousekAffineOwen');
          xx = net(sob,n);
          bds = obj.calBound;
          dx = bds(:,2) - bds(:,1);
          xCand = repmat(bds(:,1)',n,1) + xx.*repmat(dx',n,1);
      end
      
      function xCand = makeUniformsample(obj,n)
          % XCAND = MAKEUNIFORMSAMPLE(OBJ,N) returns a n-by-nVar sample matrix. The method
          % generates the sample with rand function.
          H = obj.calBound;
          nVar = size(H,1);
          dx = diff(H,[],2);
          xCand = repmat(H(:,1)',n,1)+repmat(dx',n,1).*rand(n,nVar);
      end
      
      function nVar = length(obj)
         %   NVAR = LENGTH(OBJ) returns the number of Variables in the 
         %   VariableList OBJ. 
          nVar = obj.Length;
      end
      
      function nVar = numel(obj)
          %   NVAR = NUMEL(OBJ) returns the number of Variables in the
          %   VariableList OBJ.
          nVar = obj.Length;
      end
      
      function y = checkFeasibility(obj)
          % check the feasibility of all constraints the current
          % obj VariableList have
          y = true;
          if ~isempty(obj.ExtraLinConstraint.A)
              warning('off','all');
              y = false;
              for i = 1:10
                  opt1 = optimoptions('linprog');
                  opt1.Display = 'none';
                  A = obj.ExtraLinConstraint.A;
                  ub = obj.ExtraLinConstraint.UB;
                  lb = obj.ExtraLinConstraint.LB;
                  H = obj.calBound;
                  x0 = linprog(zeros(size(A,2),1),[A;-A],[ub;-lb],[],[],H(:,1),...
                      H(:,2),opt1);
                  if ~isempty(x0)
                      y = true;
                      break
                  end
              end
          end
      end
      
      J = getJacobbian(obj,J0);
   end
end