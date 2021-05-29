classdef Dataset < handle
   % B2BDC Dataset Object
   % A dataset is a collection of models, observations, variables and their
   % respective bounds.
   %  Created: June 19, 2015   Wenyu Li
   %  Modified: June 22, 2015    Wenyu Li
   %  Modified: July 2, 2015    Wenyu Li  (option added for consistency
   %                                       measure)
   %  Modified: July 5, 2015    Wenyu Li  (Sensitivity and Range property
   %                                       for QOI added)
   %  Modified: Oct 22, 2015    Wenyu Li  (feasible point stored)
   
   properties (SetAccess = public, GetAccess = public, Hidden = true)
      FeasibleFlag = false; % Whether the feasibility includes fitting error
   end
   
   properties
      Name = '';           % Descriptive name for the dataset
      DatasetUnits = [];   % B2BDC.B2Bdataset.DatasetUnitList object
      Variables = [];      % B2BDC.B2Bvariables.VariableList object
      FeasiblePoint = [];  % Feasible point if the dataset is consistent
   end
   
   properties (Dependent)
      Length     % Number of dataset units in the dataset
   end
   
   properties (Dependent, Hidden = false)
      VarNames   % Cell array containing all variable names in the dataset
      QOINames   % Cell array containing all QOI names in the dataset
   end
   
   properties (SetAccess = private, GetAccess = public)
      ConsistencyMeasure = [];      % Lower and upper bounds to the consistency measure, [lowerBound upperBound]
      ConsistencySensitivity = [];  % Structure array of sensitivity to consistency measure
      ModelDiscrepancy = []; % Information of model discrepancy correction
   end
   
   properties (SetAccess = private, GetAccess = public, Hidden = true)
      ModelDiscrepancyFlag = false; % Whether model discrepancy is considered in the analysis
   end
   
   methods
      function obj = Dataset(dsName)
         %   OBJ = DATASET(DSNAME) returns a B2BDC.B2Bdataset.Dataset
         %   object, where DSNAME is a descriptive name for the dataset.
         if nargin > 0
            if ischar(dsName) || isstring(dsName)
               obj.Name = dsName;
            else
               error('Dataset name needs to be a character string')
            end
         end
         obj.Variables = B2BDC.B2Bvariables.VariableList();
      end
      
      function addDSunit(obj,dsUnitObj)
         % addDSunit(obj, dsUnitObj) adds a dataset unit, specified by the input dsUnitObj, to
         % the dataset
         if isempty(obj.DatasetUnits)
            obj.DatasetUnits = B2BDC.B2Bdataset.DatasetUnitList();
         end
         if isa(dsUnitObj,'B2BDC.B2Bdataset.DatasetUnit')
            for i = 1:length(dsUnitObj)
               if ~isempty(obj)
                  allName = {obj.DatasetUnits.Values.Name};
                  [~,id] = intersect(allName, dsUnitObj(i).Name);
                  if ~isempty(id)
                     error(['The dataset unit ' allName{id} ' is already exsit in the dataset'])
                  end
                  sold = obj.DatasetUnits.Values(1).ScenarioParameter;
                  if ~isempty(sold)
                     snew = dsUnitObj(i).ScenarioParameter;
                     n1 = sold.Name;
                     n2 = snew.Name;
                     [~,~,id] = intersect(n1,n2,'stable');
                     if length(id) ~= length(n1)
                        error('Dataset units have different sets of scenario parameter!');
                     elseif ~isempty(id)
                        dsUnitObj(i).ScenarioParameter.Value = dsUnitObj(i).ScenarioParameter.Value(id);
                        dsUnitObj(i).ScenarioParameter.Name = dsUnitObj(i).ScenarioParameter.Name(id);
                     end
                  end
               end
               obj.DatasetUnits = add(obj.DatasetUnits,dsUnitObj(i));
               obj.Variables = obj.Variables.addList(dsUnitObj(i).VariableList);
            end
            obj.clearConsis;
            obj.clearModelDiscrepancy;
         else
            error(['A DatasetUnit object is required as an input. Use '...
               'generateDSunit to create a DatasetUnit object before adding'])
         end
      end
      
      function y = isempty(obj)
         y = isempty(obj.DatasetUnits);
      end
      
      function y = get.Length(obj)
         y = obj.DatasetUnits.Length;
      end
      
      function y = get.VarNames(obj)
         datasetVar = obj.Variables.Values;
         n_variable = length(datasetVar);
         y = cell(n_variable,1);
         for i = 1:n_variable
            variable = datasetVar(i);
            y{i} = variable.Name;
         end
      end
      
      function y = get.QOINames(obj)
         y = {obj.DatasetUnits.Values.Name}';
      end
      
      function y = length(obj)
         y = obj.Length;
      end
      
      function y = isConsistent(obj,opt)
         % y = isConsistent(obj, opt) calculates the scalar consistency measure of the dataset 
         % with user-specified option opt, and returns dataset consistency in y
         
         % y = isConsistent(obj) calculates the scalar consistency measure of the dataset with 
         % default option, and returns dataset consistency in y
         if nargin > 1
            if ~isa(opt,'B2BDC.Option.Option')
               error(['Option input must be a B2BDC.Option object. ' ...
                  'Use generateOpt to create a B2BDC.Option object.'])
            end
         else
            opt = generateOpt;
         end
         flag1 = opt.ConsistencyMeasure;
         switch flag1
            case 'relative'
               if isempty(obj.ConsistencyMeasure)
                  if polytest(obj)
                     obj.polyConsisrel(opt);
                  elseif networktest(obj)
                     obj.networkConsisrel(opt);
                  else
                     obj.evalConsistencyrel(opt);
                  end
               end
            case 'absolute'
               if isempty(obj.ConsistencyMeasure)
                  if polytest(obj)
                     obj.polyConsisabs(opt);
                  elseif networktest(obj)
                     obj.networkConsisabs(opt);
                  else
                     obj.evalConsistencyabs(opt);
                  end
               end
         end
         if obj.ConsistencyMeasure(1) >= 0
            y = true;
         else
            y = false;
         end
         if opt.Display
            if y
               disp('The dataset is consistent')
            elseif obj.ConsistencyMeasure(2) <= 0
               disp('The dataset is inconsistent')
            else
               disp('The dataset consistency is undetermined')
            end
         end
      end
      
      function deletedUnits = deleteUnit(obj,id)
         % deletedUnits = deleteUnit(obj, id) removes the QOIs specified by the input id from
         % the dataset and retures these QOIs in deletedUnits
         deletedUnits = [];
         idx = [];
         if ischar(id) || isstring(id)
            id = {char(id)};
         end
         if iscell(id)
            allName = {obj.DatasetUnits.Values.Name};
            [~,idx] = intersect(allName, id);
         elseif isa(id, 'B2BDC.B2Bdataset.DatasetUnit')
            allName = {obj.DatasetUnits.Values.Name};
            idNames = {id.Name};
            [~,idx] = intersect(allName, idNames);
         elseif isa(id, 'double') && all(id <= obj.DatasetUnits.Length) && all(id > 0)
            idx = id;
         end
         if isempty(idx)
            error(['The ID to be removed was not found ', ...
               'or exceeded the dataset length.'])
         end
         units = obj.DatasetUnits.Values;
         obj.clearDataset;
         if ~isempty(idx)
            deletedUnits = units(idx);
            units(idx) = [];
         end
         obj.addDSunit(units);
      end
      
      function changeBound(obj,newBD,idx)
         % changeBound(obj, newBD) updates the QOI uncertainty bounds by the input
         % newBD. If newBD has 2 columns, they correspond to the lower and upper uncertainty
         % bounds; if newBD has 3 columns, the third column corresponds to the nominal value.
         %
         % changeBound(obj, newBD, idx) updates the uncertainty bounds of the QOIs
         % specified by input idx. If newBD has 2 columns, they correspond to the lower and
         % upper uncertainty bounds; if newBD has 3 columns, the third column corresponds to
         % the nominal value.
         allUnitName = {obj.DatasetUnits.Values.Name};
         if nargin == 3
            if length(idx) ~= size(newBD,1)
               error('Inconsistent new bounds size with number of dataset units')
            elseif iscell(idx)
               for i = 1:length(idx)
                  [~,id] = intersect(allUnitName, idx{i});
                  if isempty(id)
                     error(['The dataset unit ' idx{i} ' is not in the dataset'])
                  end
                  oldUnit = obj.DatasetUnits.Values(id);
                  dsName = obj.DatasetUnits.Values(id).Name;
                  newUnit = oldUnit.changeBounds(newBD(i,:));
                  obj.DatasetUnits = obj.DatasetUnits.replace('Name',dsName,newUnit);
               end
               obj.clearConsis;
            elseif isnumeric(idx)
               for i = 1:length(idx)
                  oldUnit = obj.DatasetUnits.Values(idx(i));
                  dsName = obj.DatasetUnits.Values(idx(i)).Name;
                  newUnit = oldUnit.changeBounds(newBD(i,:));
                  obj.DatasetUnits = obj.DatasetUnits.replace('Name',dsName,newUnit);
               end
               obj.clearConsis;
            else
               error('Wrong input index')
            end
         elseif nargin == 2 && size(newBD,1) == obj.Length
            for i = 1:obj.Length
               oldUnit = obj.DatasetUnits.Values(i);
               dsName = obj.DatasetUnits.Values(i).Name;
               newUnit = oldUnit.changeBounds(newBD(i,:));
               obj.DatasetUnits = obj.DatasetUnits.replace('Name',dsName,newUnit);
            end
            obj.clearConsis;
         else
            error('Wrong number of input arguments')
         end
      end
      
      function changeVarBound(obj,newBD,idx)
         % changeVarBound(obj, newBD) updates the parameter uncertainty bounds by the input
         % newBD. If newBD has 2 columns, they correspond to the lower and upper uncertainty
         % bounds; if newBD has 3 columns, the third column corresponds to the nominal value.
         %
         % changeVarBound(obj, newBD, idx) updates the uncertainty bounds of the parameters
         % specified by input idx. If newBD has 2 columns, they correspond to the lower and
         % upper uncertainty bounds; if newBD has 3 columns, the third column corresponds to
         % the nominal value.
         v1 = obj.Variables;
         v2 = v1.changeBound(newBD,idx);
         obj.Variables = v2;
      end
      
      function [units,unitIdx] = findDSunit(obj,unitName)
         % [units, unitIdx] = findDSunit(obj, unitName) finds the target QOIs specified by unitName in the dataset. The dataset units are returned in units and their indices returned
         % in unitIdx
         allName = {obj.DatasetUnits.Values.Name};
         unitIdx = find(strcmp(allName,unitName),1);
         if ~isempty(unitIdx)
            units = obj.DatasetUnits.Values(unitIdx);
         else
            units = [];
         end
      end
      
      function dsNew = clone(obj)
         % dsNew = clone(obj) makes a copy of the dataset object and returns it as dsNew.
         dsName = obj.Name;
         dsNew = B2BDC.B2Bdataset.Dataset(dsName);
         units = obj.DatasetUnits.Values;
         n = obj.Length;
         for i = 1:n
            dsNew.addDSunit(units(i));
         end
         dsNew.FeasibleFlag = obj.FeasibleFlag;
         dsNew.FeasiblePoint = obj.FeasiblePoint;
         dsNew.ConsistencyMeasure = obj.ConsistencyMeasure;
         dsNew.ConsistencySensitivity = obj.ConsistencySensitivity;
      end
      
      function bounds = calBound(obj)
         % bounds = calBound(obj) calculates the QOI uncertainty bounds.
         n = obj.Length;
         bounds = zeros(n,2);
         bounds(:,1) = [obj.DatasetUnits.Values.LowerBound]';
         bounds(:,2) = [obj.DatasetUnits.Values.UpperBound]';
      end
      
      function observes = calObserve(obj)
         % observes = calObserve(obj) calculates the measured QOI values.
         observes = [obj.DatasetUnits.Values.ObservedValue]';
      end
      
      function y = eval(obj, X, DSIdx)
         % y = eval(obj, X) evaluates all QOI outputs at input design points specified by X. Each
         % row of X corresponds to a design point. Each column of y corresponds to a QOI.
         %
         % y = eval(obj, X, DSIdx) evaluates some QOI outputs, specified by DSIdx, at input
         % design points specified by X. Each row of X corresponds to a design point. Each
         % column of y corresponds to a QOI. 
         if nargin > 0
            if nargin < 3
               DSIdx = 1:obj.Length;
            end
            if obj.ModelDiscrepancyFlag
               nMD = obj.ModelDiscrepancy.Variables.Length;
            else
               nMD = 0;
            end
            if size(X,2) ~= obj.Variables.Length && size(X,2) ~= obj.Variables.Length+nMD
               X = X';
            end
            if size(X,2) ~= obj.Variables.Length+nMD
               error('The input number of variables does not match the number of variables in the dataset')
            elseif nMD > 0
               y = obj.eval_with_discrepancy(X,DSIdx);
               return;
            end
            nSample = size(X,1);
            vObj = obj.Variables;
            if nargin == 2
               nUnit = obj.Length;
               y = zeros(nSample,nUnit);
               for i = 1:nUnit
                  tmpModel = obj.DatasetUnits.Values(i).SurrogateModel;
                  y(:,i) = tmpModel.eval(X, vObj);
               end
            elseif nargin == 3
               if isnumeric(DSIdx)
                  nUnit = length(DSIdx);
                  y = zeros(nSample, nUnit);
                  for i = 1:nUnit
                     tmpModel = obj.DatasetUnits.Values(DSIdx(i)).SurrogateModel;
                     y(:,i) = tmpModel.eval(X, vObj);
                  end
               elseif iscell(DSIdx)
                  nUnit = length(DSIdx);
                  allDS = {obj.DatasetUnits.Values.Name}';
                  y = zeros(nSample, nUnit);
                  for i = 1:nUnit
                     [~,dsID] = intersect(allDS, DSIdx{i});
                     if isempty(dsID)
                        error(['The dataset unit ' DSIdx{i} ' is not in the dataset'])
                     else
                        tmpModel = obj.DatasetUnits.Values(dsID).SurrogateModel;
                        y(:,i) = tmpModel.eval(X, vObj);
                     end
                  end
               else
                  error('Wrong input index type')
               end
            else
               error('Wrong number of input arguments')
            end
         end
      end
      
      function set.FeasiblePoint(obj,x0)
         if isempty(x0)
            obj.FeasiblePoint = [];
         else
            obj.FeasiblePoint = obj.findFeasiblePoint(x0);
         end
      end
      
      function clearConsis(obj)
         % clearConsis(obj) clears consistency-related properties of the dataset.
         obj.ConsistencyMeasure = [];
         obj.ConsistencySensitivity = [];
         obj.clearFeasiblePoint;
         obj.FeasibleFlag = false;
      end
      
      function clearModelDiscrepancy(obj)
         % clearModelDiscrepancy(obj) clears discrepancy-related properties of the dataset.
         obj.ModelDiscrepancy = [];
         obj.ModelDiscrepancyFlag = false;
         obj.clearConsis;
      end
      
      function [sv,sname] = getScenarioParameter(obj)
         % sv, sname = getScenarioParameter(obj) calculates the scenario-parameter values and
         % names of the dataset.
         ss = [obj.DatasetUnits.Values.ScenarioParameter]';
         sname = ss(1).Name;
         sv = zeros(obj.Length,length(sname));
         for i = 1:obj.Length
            sv(i,:) = ss(i).Value;
         end
      end
      
      function clearFeasiblePoint(obj)
         % clearFeasiblePoint(obj) clears the saved feasible point.
         obj.FeasiblePoint = [];
         if obj.ModelDiscrepancyFlag
            obj.ModelDiscrepancy.FeasiblePoint = [];
         end
      end
      
      function makeSubset(obj,idx)
         % makeSubset(obj, idx) updates the dataset by keeping only the QOIs specified by input
         % idx.
         if iscell(idx)
            name1 = {obj.DatasetUnits.Values.Name}';
            [~,~,idx] = intersect(idx,name1,'stable');
         end
         units = obj.DatasetUnits.Values;
         obj.clearDataset;
         if ~isempty(idx)
            remainUnits = units(idx);
         end
         obj.addDSunit(remainUnits);
      end  
   end
   
   methods (Static, Hidden = true)
      xval = q2sample(Mgd,idx,H,Xvals,V);
   end
   
   methods (Static)
      [Qall,idx] = findConicHull(Q1,Q2)
      [Qall,idx] = approxConicHull(Q1,Q2,n)
   end
   
   methods (Hidden = true)
      function clearDataset(obj)
         %   CLEARDATASET(OBJ) returns a dataset OBJ whose dataset units
         %   and variable lists have been removed.
         obj.DatasetUnits = B2BDC.B2Bdataset.DatasetUnitList();
         obj.Variables = B2BDC.B2Bvariables.VariableList();
         obj.ConsistencyMeasure = [];
         obj.ConsistencySensitivity = [];
      end
      
      function y = polytest(obj)
         %   Y = POLYTEST(OBJ) returns a logical output true if all surrogate models are
         %   polynomial models and at least one of them is not quadratic model.
         y = false;
         dsUnits = obj.DatasetUnits.Values;
         for i = 1:length(dsUnits)
            testModel = dsUnits(i).SurrogateModel;
            if isa(testModel,'B2BDC.B2Bmodels.PolyModel')
               y = true;
               continue
            elseif ~isa(testModel,'B2BDC.B2Bmodels.QModel')
               y = false;
               break;
            end
         end
      end
      
      function xx0 = findFeasiblePoint(obj,x0)
         if size(x0,1) == 1
            x0 = x0';
         end
         nVar = obj.Variables.Length;
         if obj.ModelDiscrepancyFlag
            nMD = obj.ModelDiscrepancy.Variables.Length;
            xMD = obj.ModelDiscrepancy.FeasiblePoint;
            x0 = [x0; xMD];
         else
            nMD = 0;
         end
         nx = length(x0);
         if nx == nVar+nMD
            if obj.isFeasiblePoint(x0')
               xx0 = x0(1:nVar);
               if nMD > 0
                  obj.ModelDiscrepancy.FeasiblePoint = x0(nVar+1:nVar+nMD);
               end
            else
               error('The input point is infeasible')
            end
         elseif nx == nVar
            [~,xnew] = obj.isFeasiblePoint(x0');
            if ~isempty(xnew)
               xx0 = x0;
               if nMD > 0
                  obj.ModelDiscrepancy.FeasiblePoint = xnew(nVar+1:nVar+nMD)';
               end
            else
               error('The input point is infeasible')
            end
         else
            error('The input point has a wrong dimension')
         end
      end
      
      y = eval_with_discrepancy(obj,X,DSIdx)
      evalConsistencyabs(obj,b2bopt)
      evalConsistencyDClab(obj)
      evalConsistencyrel(obj,b2bopt)
      [Qunits, Qx, Qextra, n_extra, extraIdx, L, idRQ, LBD] = getInequalQuad(obj,bds,frac)
      J = getJacobian(obj)
      directionSearch(obj,theta,x0,B2Bopt)
      [d,xopt] = calculateDistanceCurve(obj,x,opt,C)
      
      %CVX Functions
      [y,s] = cvxconsisabs(obj,yin,frac,abE)
      [y,s] = cvxconsisquadrel(obj,b2bopt,frac)
      %        [yout,sensitivity] = obj.sedumiconsisquadrel_old(b2bopt, abE);
      [minout,minSensitivity,xs] = cvxminouterbound(obj,QOIobj,frac,abE,rflag)
      [maxout,maxSensitivity,xs] = cvxmaxouterbound(obj,QOIobj,frac,abE,rflag)
      [y,s] = cvxconsisrel(obj,yin,frac,abE)
      [y,s] = cvxconsisquadabs(obj,b2bopt,abE)
      
      % Sedumi Functions
      [y,s] = sedumiconsisabs(obj,opt,abE)
      [y,s] = sedumiconsisquadabs(obj,opt,abE)
      [y,s] = sedumiconsisrel(obj,yin,opt,abE)
      [y,s] = sedumiconsisquadrel(obj,opt,abE)
      [minout,minSensitivity] = sedumiminouterbound(obj,QOIobj,frac,abE,rflag)
      [maxout,maxSensitivity] = sedumimaxouterbound(obj,QOIobj,frac,abE,rflag)
      
      % nonlinear functions
      [Qmin, Qmax, ss, xOpt, alpha] = preQOIfmincon_minB(obj,q,disflag,rflag,b2bopt,alpha)
      [Qmin, Qmax, s, xOpt, abE] = preQOIopti(obj,QOIobj,disflag,rflag,b2bopt)
      [Qmin, Qmax, s, xOpt, abE] = preQOIfmincon(obj,q,disflag,rflag,b2bopt)
      [yin_result,s,xopt,abE,flag] = relCMfmincon(obj,disflag,b2bopt)
      [yin_result,s,xopt] = relCMfminconNN(obj,disflag,b2bopt)
      [yin_result,s,xopt,abE,flag] = relCMopti(obj,disflag,b2bopt)
      [yin_result,s,xopt,abE,flag] = absCMfmincon(obj,disflag,b2bopt)
      [yin_result,s,xopt] = absCMfminconNN(obj,disflag,b2bopt)
      [yin_result,s,xopt,abE,flag] = absCMopti(obj,disflag,b2bopt)
      [vcReport, EXITFLAG] = vectorConsistencyNN(dsObj,wY,wX, nInit, opt)
      
      % samplers
      [xSample,status] = collectHRsamples_CW(obj,N,x0,opt);
      [xSample,status] = collectHRsamples(obj,N,x0,opt);
      [xSample,status] = collectAMsamples(obj,N,x0,opt);
      [xSample,status] = collectCWAMsample(obj,n,x0,opt);
      % optimization
      
      
      % Matlab Functions
      y = addlistener(obj)
      y = delete(obj)
      y = findobj(obj)
      y = findprop(obj)
      y = notify(obj)
      y = le(obj)
      y = lt(obj)
      y = ne(obj)
      y = ge(obj)
      y = gt(obj)
      y = eq(obj)
      
      
   end
   
end
