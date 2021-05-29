function [vcReport, EXITFLAG] = vectorConsistencyNN(dsObj,wY,wX, opt)
% vectorConsistency   Evaluates the vector consistency of a given dataset
% with given weights
%    V = obj.vectorConsistency(wY,wX)  returns a structure containing the
%    results of the vector consistency analysis, including the necessary
%    relaxations to the QOI bounds and the variable bounds to render the
%    dataset consistent.
%
% The inputs are:
% --------------------------------------------------------------------------
%     dsObj - An inconsistent dataset
%     wY    - A obj.Length-by-2 matrix of weights, where the first column
%            acts on the relaxations to the corresponding QOI lower bounds
%            and the second column acts on the QOI upper bounds. A weight
%            of zero enforces corresponding QOI constraint, i.e. allows no
%            relaxation.
%     wX    - A obj.Variables.Length-by-2 matrix of weights to the
%            relaxations on the
%            variables.
%     nInit - number of initial conditions for fmincon, sampled from SDP
%            result. In some cases, multiple initial conditions will
%            improve the quality of the returned solution.
%     opt   - a B2BDC option
% --------------------------------------------------------------------------
% Instead of specifying a weight matrix, setting wY and/or wX to the
% character strings 'perc','uwidth','unit', or 'null' automatically activates one
% of the standard weighting configurations (percentage change, uncertainty
% width, absolute change, or zero)


%Step 1: Process inputs
nVar = dsObj.Variables.Length; %number of variables
nExp = dsObj.DatasetUnits.Length; %number of experiments
yBnds = dsObj.calBound; %[lb ub measurement];
xBnds = dsObj.Variables.calBound;

if ischar(wY) && strcmpi(wY, 'perc') %percentage consistency measure
   wY = abs(yBnds);
elseif ischar(wY) && strcmpi(wY, 'uwidth') %symmetric uncertainy width
   wY = repmat((yBnds(:,2) - yBnds(:,1)),1,2);
elseif ischar(wY) && strcmpi(wY, 'unit')
   wY = ones(nExp,2);
elseif ischar(wY) && strcmpi(wY, 'null') || norm(wY) == 0
   wY = zeros(nExp,2);
elseif size(wY,1) ~= nExp
   error('Improper QOI weight specification')
end
if ischar(wX) && strcmpi(wX, 'perc')
   wX = abs(xBnds);
elseif ischar(wX) && strcmpi(wX, 'uwidth') %symmetric uncertainty width
   wX = repmat((xBnds(:,2) - xBnds(:,1)),1,2);
elseif ischar(wX) && strcmpi(wX, 'unit')
   wX = ones(nVar,2);
elseif ischar(wX) && strcmpi(wX, 'null') || norm(wX)==0
   wX = zeros(nVar,2);
elseif size(wX,1) ~= nVar
   error('Improper variable weight specification')
end

if nargin < 4
   opt = generateOpt;
end
nInit = opt.OptimOption.RandomStart;
units = dsObj.DatasetUnits.Values;
absErr = zeros(nExp,1);
if opt.AddFitError
   for j = 1:nExp
      if ~isempty(units(j).SurrogateModel.ErrorStats)
         absErr(j) = units(j).SurrogateModel.ErrorStats.absMax;
      end
   end
end

fminopt = optimoptions('fmincon','Display','none','GradObj','on',...
   'GradConstr','on','Algorithm','interior-point','MaxIter',5000,...
   'TolFun',1e-6,'TolCon',1e-6,'DerivativeCheck','off','MaxFunctionEvaluations',1e4,...
   'FiniteDifferenceType','central');

% Ignore noninfluential/redundant relaxations, i.e. those with zero weight.
% Used to build the problem with the minimum number of dimensions needed.

indYLb = find(wY(:,1)>0);
indYUb = find(wY(:,2)>0);

nYLb = numel(indYLb);
nYUb = numel(indYUb);

nVCM = nYLb+nYUb;

wyl = wY(indYLb,1);
wyu = wY(indYUb,2);

vList = dsObj.Variables;
if ~isempty(vList.ExtraLinConstraint.A)
   tolerance = 1e-5;
   A0 = vList.ExtraLinConstraint.A;
   n3 = size(A0,1);
   Aub = vList.ExtraLinConstraint.UB;
   Alb = vList.ExtraLinConstraint.LB;
   eqTest = (Aub-Alb)./sum(abs(A0),2);
   dA = Aub - Alb;
   ieq = (eqTest <= tolerance);
   if any(ieq)
      A1 = [A0(~ieq,:); -A0(~ieq,:)];
      B1 = [Aub(~ieq); -Alb(~ieq)];
      Aeq = A0(ieq,:);
      Beq = 0.5*(Alb(ieq)+Aub(ieq));
   else
      A1 = [A0; -A0];
      Aeq = [];
      B1 = [Aub; -Alb];
      Beq = [];
   end
   n1 = size(A1,1);
   n2 = size(Aeq,1);
   A1 = [zeros(n1,1), A1, zeros(n1,nVCM)];
   if n2 ~= 0
      Aeq = [zeros(n2,1), Aeq, zeros(nVCM)];
   end
else
   A1 = [];
   B1 = [];
   Aeq = [];
   Beq = [];
end
allVarnames = dsObj.VarNames;
idall = cell(nExp,1);
fg = cell(nExp,1);
for j = 1:nExp
   tmodel = units(j).SurrogateModel;
   [~,~,idall{j}] = intersect(tmodel.VarNames,allVarnames,'stable');
   fg{j} = tmodel.NetGradientFcn;
end
if opt.Display
   disp('=======================================================');
   disp('Searching inner bound...');
   disp('=======================================================');
end
xStart = [vList.makeLHSsample(nInit) zeros(nInit,nVCM)];
xOpt = cell(nInit,1);
yOpt = zeros(nInit,1);
flags = zeros(nInit,1);
for i = 1:nInit
   [xOpt{i},yOpt(i),flags(i)] = fmincon(@funmin,[0;xStart(i,:)'],...
      A1,B1,Aeq,Beq,[0;xBnds(:,1);zeros(nVCM,1)],[inf;xBnds(:,2);inf(nVCM,1)],@neq,fminopt);
end
xOpt(flags<=0) = [];
yOpt(flags<=0) = [];

if isempty(yOpt)
   vcReport = [];
   EXITFLAG = -1;
else
   [~,ID] = min(yOpt);
   zFeas = xOpt{ID}(2:end);
   DeltaLb = zeros(nExp,1);
   DeltaLb(indYLb) = zFeas(nVar+1:nVar+nYLb);
   DeltaUb = zeros(nExp,1);
   DeltaUb(indYUb) = zFeas(nVar+nYLb+1:nVar+nYLb+nYUb);
   deltaLb = zeros(nVar,1);
   deltaUb = zeros(nVar,1);
   vcReport.Objective.SDPLowerBound = [];
   vcReport.Objective.Feasible = xOpt{ID}(1);
   vcReport.FeasiblePoint = zFeas(1:nVar);
   vcReport.Relaxations.yLowerBound = DeltaLb;
   vcReport.Relaxations.yUpperBound = DeltaUb;
   vcReport.Relaxations.xLowerBound = deltaLb;
   vcReport.Relaxations.xUpperBound = deltaUb;
   vcReport.Weights.yLowerBound = wY(:,1);
   vcReport.Weights.yUpperBound = wY(:,2);
   vcReport.Weights.xLowerBound = zeros(nVar,1);
   vcReport.Weights.xUpperBound = zeros(nVar,1);
   EXITFLAG = 1;
end

   function [y,gy] = funmin(x)
      y = x(1);
      gy = zeros(nVar+nVCM+1,1);
      gy(1) = 1;
   end

   function [c,ceq,g,geq] = neq(x)
      c = zeros(2*nExp+1,1);
      g = zeros(nVar+nVCM+1,2*nExp+1);
      for ii = 1:nExp
         id = idall{ii}+1;
         c(2*ii-1) = units(ii).SurrogateModel.eval(x(id)');
         g(id,2*ii-1) = fg{ii}(x(id));
      end
      c(2:2:end-1) = -c(1:2:end-1);
      g(:,2:2:end-1) = -g(:,1:2:end-1);
      c(1:2:end-1) = c(1:2:end-1)-yBnds(:,2);
      c(2:2:end-1) = yBnds(:,1)+c(2:2:end-1);
      c(end) = sum(x(nVar+2:end))-x(1);
      g(1,end) = -1;
      g(nVar+2:end,end) = 1;
      c(2*indYLb) = c(2*indYLb)-wyl.*x(nVar+2:nVar+nYLb+1);
      g(nVar+2:nVar+nYLb+1,2*indYLb) = -diag(wyl);
      c(2*indYUb-1) = c(2*indYUb-1)-wyu.*x(nVar+nYLb+2:nVar+nYLb+nYUb+1);
      g(nVar+nYLb+2:nVar+nYLb+nYUb+1,2*indYUb-1) = -diag(wyu);
      ceq = [];
      geq = [];
   end
end


