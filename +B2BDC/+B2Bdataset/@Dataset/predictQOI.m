function [QOIrange, QOISensitivity, xOpt] = predictQOI(obj, QOIobj, B2Bopt, QOIcorrection)
% Returns the inner and outer bounds, the feasible points of maximum and minimum 
% values of the QOI model subject to the constraints that both dataset units and
% variables are within their respective bounds. The sensitivity of the
% outer bounds with respect to experimental bounds and variable bounds is
% also returned.
% The input arguments are:
%    QOIobj - A B2BDC.B2Bmodel.Model object that defines the QOI
%    B2Bopt - B2BDC.Option object (optional)
%    QOIcorrection - structure with GroupIndex and Value (scenario value)
% The output are:
%    QOIrange - A structure defines the range of QOI subject to the dataset
%    Xopt - A nVar-by-2 double matrix
%
% [QOIrange, QOISensitivity, xOpt] = predictQOI(obj, QOIobj, B2Bopt) predicts the
% uncertainty interval of the QOI specified by QOIobj using the dataset obj. The option
% variable is specified by B2Bopt. The calculated uncertainty interval, sensitivity, and
% optimal points for the inner bounds, are saved in QOIrange, QOISensitivity, and xOpt,
% respectively.
%
% [QOIrange, QOISensitivity, xOpt] = predictQOI(obj, QOIobj, B2Bopt, QOIcorrection)
% predicts the uncertainty interval of the QOI specified by QOIobj using the dataset obj.
% Model discrepancy is considered in the analysis and the corresponding information of
% the QOI is specified by QOIcorrection. The option variable is specified by B2Bopt. The
% calculated uncertainty interval, sensitivity, and optimal points for the inner bounds,
% are saved in QOIrange, QOISensitivity and xOpt, respectively

% Created: June 15, 2015   Wenyu Li
%  Modified: July 5, 2015   Wenyu Li (Sensitivity added)
%  Modified: July 14, 2015   Wenyu Li  (Gradient and hessian added)

if nargin < 3
   B2Bopt = generateOpt;
elseif ~isa(B2Bopt,'B2BDC.Option.Option')
   error('Wrong option object')
end
if nargin > 3 && ~isempty(QOIcorrection)
   q.Model = QOIobj;
   q.Correction = QOIcorrection;
else
   q.Model = QOIobj;
   q.Correction.GroupIndex = 0;
   q.Correction.Value = [];
end
name1 = QOIobj.VarNames;
name2 = obj.VarNames;
if obj.ModelDiscrepancyFlag
   name2 = [name2; {obj.ModelDiscrepancy.Variables.Values.Name}'];
end
[~,id] = intersect(name1,name2);
if length(id) < length(name1)
   rflag = true;
else
   rflag = false;
end
if ~obj.isConsistent(B2Bopt)
   error('The dataset is inconsistent')
else
   tmpCM = obj.ConsistencyMeasure;
   tmpS = obj.ConsistencySensitivity;
   tmpFlag = obj.FeasibleFlag;
   tmpx0 = obj.FeasiblePoint;
end
% xOpt = zeros(obj.Variables.Length,2);
pFlag = B2Bopt.Prediction;
if rflag
   tVar = QOIobj.Variables;
   ntVar = tVar.Length;
   tCoef = zeros(ntVar+1);
   tCoef(1,1) = 1;
   tM = B2BDC.B2Bmodels.QModel(tCoef,tVar);
   tUnit = generateDSunit('redundant',tM,[0,2]);
   obj.addDSunit(tUnit);
   if ~obj.isConsistent(B2Bopt)
      error('The dataset is inconsistent')
   end
end
polyflag = polytest(obj);
nnflag = networktest(obj);
if polyflag
   if B2Bopt.Display
      disp('=======================================================');
      disp('Calculating QOI range...');
      disp('=======================================================');
   end
   [QOIrange, QOISensitivity, xOpt] = obj.preQOIpoly(QOIobj,B2Bopt,rflag);
   if B2Bopt.Display
      disp('The calculation is done')
      disp(['Minimum value of QOI is within: [' num2str(QOIrange.min(1)) ' ' num2str(QOIrange.min(2)) ']'])
      disp(['Maximum value of QOI is within: [' num2str(QOIrange.max(1)) ' ' num2str(QOIrange.max(2)) ']'])
   end
   return
elseif nnflag
   if B2Bopt.Display
      disp('=======================================================');
      disp('Calculating QOI range...');
      disp('=======================================================');
   end
   [QOIrange, QOISensitivity, xOpt] = obj.preQOInn(q,B2Bopt,rflag);
   if B2Bopt.Display
      disp('The calculation is done')
      disp(['Minimum value of QOI is within: [' num2str(QOIrange.min(1)) ' ' num2str(QOIrange.min(2)) ']'])
      disp(['Maximum value of QOI is within: [' num2str(QOIrange.max(1)) ' ' num2str(QOIrange.max(2)) ']'])
   end
   return
end

if ~strcmp(pFlag,'outer')
%    tic;
   [minin, maxin, s, xOpt, abE] = obj.preQOIfmincon(q,B2Bopt.Display,rflag,B2Bopt);
%    [minin, maxin, s, xOpt, abE] = obj.preQOIopti(q,B2Bopt.Display,rflag,B2Bopt);
%    toc;
   QOISensitivity.Inner = s;
else
   units = obj.DatasetUnits.Values;
   abE = zeros(length(units),1);
   if B2Bopt.AddFitError
      for j = 1:length(units)
         if ~isempty(units(j).SurrogateModel.ErrorStats.absMax)
            abE(j) = units(j).SurrogateModel.ErrorStats.absMax;
         end
      end
   end
   QOISensitivity.Inner = [];
end

if ~strcmp(pFlag,'inner')
   if B2Bopt.Display
      disp('=======================================================');
      disp('Calculating outer bound...');
      disp('=======================================================');
   end
   frac = B2Bopt.ExtraLinFraction;
%    [minout,minSens] = obj.sedumiminouterbound(q,frac,abE,rflag);
%    [maxout,maxSens] = obj.sedumimaxouterbound(q,frac,abE,rflag);
   warning('off','all');
   [minout,minSens,xs_min] = obj.cvxminouterbound(q,frac,abE,rflag);
   [maxout,maxSens,xs_max] = obj.cvxmaxouterbound(q,frac,abE,rflag);
   warning('on','all');
   minSens.expu = -minSens.expu;
   minSens.expl = -minSens.expl;
   minSens.varu = -minSens.varu;
   minSens.varl = -minSens.varl;
   minSens.linu = -minSens.linu;
   minSens.linl = -minSens.linl;
%    maxSens.expu = maxSens.expu;
%    maxSens.expl = maxSens.expl;
%    maxSens.varu = maxSens.varu;
%    maxSens.varl = maxSens.varl;
%    maxSens.linu = maxSens.linu;
%    maxSens.linl = maxSens.linl;
   QOISensitivity.Outer.min = minSens;
   QOISensitivity.Outer.max = maxSens;
else
   QOISensitivity.Outer = [];
end
switch pFlag
   case 'both'
      QOIrange.min = [minout, minin];
      QOIrange.max = [maxin, maxout];
      if B2Bopt.Display
         disp('The calculation is done')
         disp(['Minimum value of QOI is within: [' num2str(minout) ' ' num2str(minin) ']'])
         disp(['Maximum value of QOI is within: [' num2str(maxin) ' ' num2str(maxout) ']'])
      end
   case 'inner'
      QOIrange.min = minin;
      QOIrange.max = maxin;
      if B2Bopt.Display
         disp('The calculation is done')
         disp(['Minimum value (inner bound) of QOI is ' num2str(minin)])
         disp(['Maximum value (inner bound) of QOI is ' num2str(maxin)])
      end
   case 'outer'
      QOIrange.min = minout;
      QOIrange.max = maxout;
      if B2Bopt.Display
         disp('The calculation is done')
         disp(['Minimum value (outer bound) of QOI is ' num2str(minout)])
         disp(['Maximum value (outer bound) of QOI is ' num2str(maxout)])
      end
end
if rflag
   obj.deleteUnit(obj.Length);
end

obj.ConsistencyMeasure = tmpCM;
obj.ConsistencySensitivity = tmpS;
obj.FeasibleFlag = tmpFlag;
obj.FeasiblePoint = tmpx0;

end

