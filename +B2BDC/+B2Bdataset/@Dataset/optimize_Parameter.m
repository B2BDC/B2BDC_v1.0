function [xopt,yopt] = optimize_Parameter(obj,logFlag,weight,opt)
   % [xopt, yopt] = optimize Parameter(obj, logFlag, weight) calculates an optimized parameter 
   % vector based on a specified criteria. The logFlag indicates whether the response
   % surface of the corresponding QOI is in the log scale. The minimized object function is
   % specified by the option variable and the weights asscoiated with individual terms in the
   % object function is specified by weight. In this case, the option variable is the default
   % B2BDC option. The calculated optimal parameter vector and minimized object value
   % are saved in xopt and yopt.
   %
   % [xopt, yopt] = optimize Parameter(obj, logFlag, weight, opt) calculates an optimized parameter 
   % vector based on a specified criteria. The logFlag indicates whether the response
   % surface of the corresponding QOI is in the log scale. The minimized object function is
   % specified by the option variable and the weights asscoiated with individual terms in the
   % object function is specified by weight. In this case, the option variable is specified by
   % opt. The calculated optimal parameter vector and minimized object value are saved
   % in xopt and yopt.

%  Created: Oct 25, 2017

if nargin < 4
   opt = generateOpt;
end
optimOpt = opt.OptimOption;
method = optimOpt.OptimizationMethod;
logType = strcmp(optimOpt.Logtype,'nature');
e10 = log(10);
vList = obj.Variables;
nVar = vList.Length;
varNom = [vList.Values.NominalValue]';
% varNom = zeros(nVar,1);
ns = optimOpt.RandomStart;
H = obj.Variables.calBound;
LB = H(:,1);
UB = H(:,2);
if obj.ModelDiscrepancyFlag
   MDvar = obj.ModelDiscrepancy.Variables;
   nMD = MDvar.Length;
   MDbd = MDvar.calBound;
   LB = [LB; MDbd(:,1)];
   UB = [UB; MDbd(:,2)];
   varNom = [varNom; obj.ModelDiscrepancy.Variables.calNominal];
else
   nMD = 0;
end
if isempty(logFlag)
   logFlag = false(obj.Length+nVar+nMD,1);
elseif numel(logFlag) ~= obj.Length+nVar+nMD
   error('The input logFlag has a wrong dimension');
end
if ~strcmp(method,'LSH')
   if ~obj.isConsistent(opt)
      disp('The dataset is inconsistent')
      xopt = [];
      yopt = [];
      return
   end
   xFea = obj.FeasiblePoint;
   if nMD > 0
      xFea = [xFea; obj.ModelDiscrepancy.FeasiblePoint];
   end
   step0 = opt.SampleOption.StepInterval;
   opt.SampleOption.StepInterval = max(step0, 1e2);
   xx = obj.collectHRsamples_CW(ns,xFea',opt);
   opt.SampleOption.StepInterval = step0;
   x0 = xx.x;
   if ~all(obj.isFeasiblePoint(x0))
      error('Some sampled points are infeasible')
   end
   if strcmp(method,'1NF')
      optimOpt.PenaltyWeight = 'user-defined';
   end
else
%    x0 = vList.makeLHSsample(ns);
   x0 = obj.Variables.makeDesignSample(ns);
   if nMD > 0
%       x0 = [x0 MDvar.makeLHSsample(ns)];
      x0 = [x0 zeros(ns,nMD)];
   end
end

% varNom = zeros(nVar,1);
if strcmp(method,'LSH')
   fminopt = optimoptions('fmincon','Display','none','GradObj','on','MaxFunctionEvaluations',5e4,...
      'GradConstr','on','Algorithm','interior-point','MaxIter',1e4,'TolFun',1e-10,'TolCon',1e-10,...
      'StepTolerance',1e-20,'DerivativeCheck','off','FiniteDifferenceType','central');
else
   fminopt = optimoptions('fmincon','Display','none','GradObj','on','MaxFunctionEvaluations',1e4,...
      'GradConstr','on','Algorithm','interior-point','MaxIter',5e3,'TolFun',1e-10,'TolCon',1e-10,...
      'StepTolerance',1e-20,'DerivativeCheck','off','FiniteDifferenceType','central');
end
if ~isempty(vList.ExtraLinConstraint.A)
   tolerance = 1e-5;
   A0 = vList.ExtraLinConstraint.A;
   Aub = vList.ExtraLinConstraint.UB;
   Alb = vList.ExtraLinConstraint.LB;
   eqTest = (Aub - Alb)./sum(abs(A0),2);
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
else
   A1 = [];
   B1 = [];
   Aeq = [];
   Beq = [];
end
bds = obj.calBound;
obs = obj.calObserve;
units = obj.DatasetUnits.Values;
n_units = length(units);
if logType
   obs(logFlag(1:n_units)) = exp(obs(logFlag(1:n_units)));
   varNom(logFlag(n_units+1:end)) = exp(varNom(logFlag(n_units+1:end)));
else
   obs(logFlag(1:n_units)) = 10.^(obs(logFlag(1:n_units)));
   varNom(logFlag(n_units+1:end)) = 10.^(varNom(logFlag(n_units+1:end)));
end
abE = zeros(length(units),1);
if opt.AddFitError
   for j = 1:length(units)
      if ~isempty(units(j).SurrogateModel.ErrorStats.absMax)
         abE(j) = units(j).SurrogateModel.ErrorStats.absMax;      
      end
   end
end
% abE = abE-0.5*diff(bds,[],2)*tolY;
switch optimOpt.PenaltyWeight
   case 'relative'
      if nargin < 4
         w = [1./obs; zeros(nVar+nMD,1)];
      else
         w = weight;
      end
   case 'absolute'
      if nargin < 4
         w = [ones(n_units,1); zeros(nVar+nMD,1)];
      else
         w = weight;
      end
   case 'user-defined'
      if ~isempty(weight)
         if length(weight) ~= n_units + nVar + nMD
            error('Dimension of input weight is wrong')
         else
            w = weight;
         end
      elseif strcmp(optimOpt.OptimizationMethod,'1NF')
         w = [zeros(n_units,1); ones(nVar,1); zeros(nMD,1)];
      else
         error('Weights are not specified');
      end
end
[idall,Qall,Nall,Dall] = obj.getQ_RQ_expansion;
if ~isempty(A1)
   A1 = [A1 zeros(size(A1,1),nMD)]; 
end
if ~isempty(Aeq)
   Aeq = [Aeq zeros(size(Aeq,1),nMD)];
end
bds = bds+[-abE abE];

xopt = zeros(nVar+nMD,ns);
yopt = zeros(1,ns);
exitflag = zeros(1,ns);
if logType
   for j = 1:ns
      if strcmp(optimOpt.PenaltyWeight,'relative')
         switch method
            case 'LSF'
               [xopt(:,j),yopt(j),exitflag(j)] = fmincon(@min_LS_rel,x0(j,:)',...
                  A1,B1,Aeq,Beq,LB,UB,@neq_F,fminopt);
            case 'LSH'
               [xopt(:,j),yopt(j),exitflag(j)] = fmincon(@min_LS_rel,x0(j,:)',...
                  A1,B1,Aeq,Beq,LB,UB,[],fminopt);
            case '1NF'
               [xopt(:,j),yopt(j),exitflag(j)] = fmincon(@min_1N,x0(j,:)',...
                  A1,B1,Aeq,Beq,LB,UB,@neq_F,fminopt);
         end
      else
         switch method
            case 'LSF'
               [xopt(:,j),yopt(j),exitflag(j)] = fmincon(@min_LS,x0(j,:)',...
                  A1,B1,Aeq,Beq,LB,UB,@neq_F,fminopt);
            case 'LSH'
               [xopt(:,j),yopt(j),exitflag(j)] = fmincon(@min_LS,x0(j,:)',...
                  A1,B1,Aeq,Beq,LB,UB,[],fminopt);
            case '1NF'
               [xopt(:,j),yopt(j),exitflag(j)] = fmincon(@min_1N,x0(j,:)',...
                  A1,B1,Aeq,Beq,LB,UB,@neq_F,fminopt);
         end
      end
   end
else
   for j = 1:ns
      if strcmp(optimOpt.PenaltyWeight,'relative')
         switch method
            case 'LSF'
               [xopt(:,j),yopt(j),exitflag(j)] = fmincon(@min_LS2_rel,x0(j,:)',...
                  A1,B1,Aeq,Beq,LB,UB,@neq_F,fminopt);
            case 'LSH'
               [xopt(:,j),yopt(j),exitflag(j)] = fmincon(@min_LS2_rel,x0(j,:)',...
                  A1,B1,Aeq,Beq,LB,UB,[],fminopt);
            case '1NF'
               [xopt(:,j),yopt(j),exitflag(j)] = fmincon(@min_1N2,x0(j,:)',...
                  A1,B1,Aeq,Beq,LB,UB,@neq_F,fminopt);
         end
      else
         switch method
            case 'LSF'
               [xopt(:,j),yopt(j),exitflag(j)] = fmincon(@min_LS2,x0(j,:)',...
                  A1,B1,Aeq,Beq,LB,UB,@neq_F,fminopt);
            case 'LSH'
               [xopt(:,j),yopt(j),exitflag(j)] = fmincon(@min_LS2,x0(j,:)',...
                  A1,B1,Aeq,Beq,LB,UB,[],fminopt);
            case '1NF'
               [xopt(:,j),yopt(j),exitflag(j)] = fmincon(@min_1N2,x0(j,:)',...
                  A1,B1,Aeq,Beq,LB,UB,@neq_F,fminopt);
         end
      end
   end
end
id = find(exitflag<=0);
xopt(:,id) = [];
yopt(id) = [];
if isempty(yopt)
   disp('Optimal point is not found.')
else
   [~,id] = min(yopt);
   xopt = xopt(:,id);
   yopt = yopt(id);
   if ~strcmp(method,'1NF')
      yopt = yopt/n_units;
   end
end

   function [y,gy] = min_LS(x)
      y = 0;
      gy = zeros(nVar+nMD,1);
      unitID = find(w(1:n_units));
      for i = unitID'
         tw = w(i);
         id1 = idall{i};
         if isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.RQModel')
            N = Nall{i};
            D = Dall{i};
            tmpN = [1;x(id1)]'*N*[1;x(id1)];
            tmpD = [1;x(id1)]'*D*[1;x(id1)];
            if logFlag(i)
               y = y+tw*(exp(tmpN/tmpD)-obs(i))^2;
            else
               y = y+tw*(tmpN/tmpD-obs(i))^2;
            end
            if logFlag(i)
               gy(id1) = gy(id1)+2*tw*(exp(tmpN/tmpD)-obs(i))*exp(tmpN/tmpD)...
                  /tmpD^2*(tmpD*(N(2:end,2:end)*x(id1)+N(2:end,1))-...
                  tmpN*(D(2:end,2:end)*x(id1)+D(2:end,1)));
            else
               gy(id1) = gy(id1)+2*tw*(tmpN/tmpD-obs(i))/tmpD^2*...
                  (tmpD*(N(2:end,2:end)*x(id1)+N(2:end,1))-...
                  tmpN*(D(2:end,2:end)*x(id1)+D(2:end,1)));
            end
         elseif isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.QModel')
            quadCoef = Qall{i};
            tmpQ = [1;x(id1)]'*quadCoef*[1;x(id1)];
            if logFlag(i)
               y = y + tw*(exp(tmpQ) - obs(i))^2;
            else
               y = y + tw*(tmpQ - obs(i))^2;
            end
            if logFlag(i)
               gy(id1) = gy(id1) + 4*tw*(exp(tmpQ)-obs(i))*exp(tmpQ)*...
                  (quadCoef(2:end,2:end)*x(id1)+quadCoef(2:end,1));
            else
               gy(id1) = gy(id1) + 4*tw*(tmpQ-obs(i))*...
                  (quadCoef(2:end,2:end)*x(id1)+quadCoef(2:end,1));
            end
         end
      end
      if strcmp(optimOpt.PenaltyWeight,'user-defined')
         varID = find(w(n_units+1:end));
         for i = varID'
            tw = w(i+n_units);
            if logFlag(i+n_units)
               y = y + tw*(exp(x(i))-varNom(i))^2;
               gy(i) = gy(i) + 2*tw*(exp(x(i))-varNom(i))*exp(x(i));
            else
               y = y + tw*(x(i)-varNom(i))^2;
               gy(i) = gy(i) + 2*tw*(x(i)-varNom(i));
            end
         end
      end
   end

   function [y,gy] = min_LS_rel(x)
      y = 0;
      gy = zeros(nVar+nMD,1);
      unitID = find(w(1:n_units));
      for i = unitID'
         tw = w(i);
         id1 = idall{i};
         if isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.RQModel')
            N = Nall{i};
            D = Dall{i};
            tmpN = [1;x(id1)]'*N*[1;x(id1)];
            tmpD = [1;x(id1)]'*D*[1;x(id1)];
            if logFlag(i)
               y = y+abs(tw*(exp(tmpN/tmpD)-obs(i)));
            else
               y = y+abs(tw*((tmpN/tmpD-obs(i))));
            end
            if logFlag(i)
               gy(id1) = gy(id1)+sign(tw*((exp(tmpN/tmpD)-obs(i))))*tw*exp(tmpN/tmpD)...
                  /tmpD^2*(tmpD*(N(2:end,2:end)*x(id1)+N(2:end,1))-...
                  tmpN*(D(2:end,2:end)*x(id1)+D(2:end,1)));
            else
               gy(id1) = gy(id1)+sign(tw*(tmpN/tmpD-obs(i)))*tw/tmpD^2*...
                  (tmpD*(N(2:end,2:end)*x(id1)+N(2:end,1))-...
                  tmpN*(D(2:end,2:end)*x(id1)+D(2:end,1)));
            end
         elseif isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.QModel')
            quadCoef = Qall{i};
            tmpQ = [1;x(id1)]'*quadCoef*[1;x(id1)];
            if logFlag(i)
               y = y + abs(tw*(exp(tmpQ) - obs(i)));
            else
               y = y + abs(tw*(tmpQ - obs(i)));
            end
            if logFlag(i)
               gy(id1) = gy(id1) + 2*sign(tw*(exp(tmpQ)-obs(i)))*tw*exp(tmpQ)*...
                  (quadCoef(2:end,2:end)*x(id1)+quadCoef(2:end,1));
            else
               gy(id1) = gy(id1) + 2*sign(tw*(tmpQ-obs(i)))*tw*...
                  (quadCoef(2:end,2:end)*x(id1)+quadCoef(2:end,1));
            end
         end
      end
      if strcmp(optimOpt.PenaltyWeight,'user-defined')
         varID = find(w(n_units+1:end));
         for i = varID'
            tw = w(i+n_units);
            if logFlag(i+n_units)
               y = y + tw*(exp(x(i))-varNom(i))^2;
               gy(i) = gy(i) + 2*tw*(exp(x(i))-varNom(i))*exp(x(i));
            else
               y = y + tw*(x(i)-varNom(i))^2;
               gy(i) = gy(i) + 2*tw*(x(i)-varNom(i));
            end
         end
      end
   end

   function [y,gy] = min_LS2(x)
      y = 0;
      gy = zeros(nVar+nMD,1);
      unitID = find(w(1:n_units));
      for i = unitID'
         tw = w(i);
         id1 = idall{i};
         if isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.RQModel')
            N = Nall{i};
            D = Dall{i};
            tmpN = [1;x(id1)]'*N*[1;x(id1)];
            tmpD = [1;x(id1)]'*D*[1;x(id1)];
            if logFlag(i)
               y = y+tw*(10^(tmpN/tmpD)-obs(i))^2;
            else
               y = y+tw*(tmpN/tmpD-obs(i))^2;
            end
            if logFlag(i)
               gy(id1) = gy(id1)+2*tw*e10*(10^(tmpN/tmpD)-obs(i))*10^(tmpN/tmpD)...
                  /tmpD^2*(tmpD*(N(2:end,2:end)*x(id1)+N(2:end,1))-...
                  tmpN*(D(2:end,2:end)*x(id1)+D(2:end,1)));
            else
               gy(id1) = gy(id1)+2*tw*(tmpN/tmpD-obs(i))/tmpD^2*...
                  (tmpD*(N(2:end,2:end)*x(id1)+N(2:end,1))-...
                  tmpN*(D(2:end,2:end)*x(id1)+D(2:end,1)));
            end
         elseif isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.QModel')
            quadCoef = Qall{i};
            tmpQ = [1;x(id1)]'*quadCoef*[1;x(id1)];
            if logFlag(i)
               y = y + tw*(10^(tmpQ) - obs(i))^2;
            else
               y = y + tw*(tmpQ - obs(i))^2;
            end
            if logFlag(i)
               gy(id1) = gy(id1) + 4*tw*e10*(10^(tmpQ)-obs(i))*10^(tmpQ)*...
                  (quadCoef(2:end,2:end)*x(id1)+quadCoef(2:end,1));
            else
               gy(id1) = gy(id1) + 4*tw*(tmpQ-obs(i))*...
                  (quadCoef(2:end,2:end)*x(id1)+quadCoef(2:end,1));
            end
         end
      end
      if strcmp(optimOpt.PenaltyWeight,'user-defined')
         varID = find(w(n_units+1:end));
         for i = varID'
            tw = w(i+n_units);
            if logFlag(i+n_units)
               y = y + tw*(10^(x(i))-varNom(i))^2;
               gy(i) = gy(i) + 2*tw*e10*(10^(x(i))-varNom(i))*10^(x(i));
            else
               y = y + tw*(x(i)-varNom(i))^2;
               gy(i) = gy(i) + 2*tw*(x(i)-varNom(i));
            end
         end
      end
   end

   function [y,gy] = min_LS2_rel(x)
      y = 0;
      gy = zeros(nVar+nMD,1);
      unitID = find(w(1:n_units));
      for i = unitID'
         tw = w(i);
         id1 = idall{i};
         if isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.RQModel')
            N = Nall{i};
            D = Dall{i};
            tmpN = [1;x(id1)]'*N*[1;x(id1)];
            tmpD = [1;x(id1)]'*D*[1;x(id1)];
            if logFlag(i)
               y = y+abs(tw*(10^(tmpN/tmpD)-obs(i)));
            else
               y = y+abs(tw*(tmpN/tmpD-obs(i)));
            end
            if logFlag(i)
               gy(id1) = gy(id1)+tw*e10*sign(tw*(10^(tmpN/tmpD)-obs(i)))*10^(tmpN/tmpD)...
                  /tmpD^2*(tmpD*(N(2:end,2:end)*x(id1)+N(2:end,1))-...
                  tmpN*(D(2:end,2:end)*x(id1)+D(2:end,1)));
            else
               gy(id1) = gy(id1)+tw*sign(tw*(tmpN/tmpD-obs(i)))/tmpD^2*...
                  (tmpD*(N(2:end,2:end)*x(id1)+N(2:end,1))-...
                  tmpN*(D(2:end,2:end)*x(id1)+D(2:end,1)));
            end
         elseif isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.QModel')
            quadCoef = Qall{i};
            tmpQ = [1;x(id1)]'*quadCoef*[1;x(id1)];
            if logFlag(i)
               y = y + abs(tw*(10^(tmpQ) - obs(i)));
            else
               y = y + abs(tw*(tmpQ - obs(i)));
            end
            if logFlag(i)
               gy(id1) = gy(id1) + 2*tw*e10*sign(tw*(10^(tmpQ)-obs(i)))*10^(tmpQ)*...
                  (quadCoef(2:end,2:end)*x(id1)+quadCoef(2:end,1));
            else
               gy(id1) = gy(id1) + 2*tw*sign(tw*(tmpQ-obs(i)))*...
                  (quadCoef(2:end,2:end)*x(id1)+quadCoef(2:end,1));
            end
         end
      end
      if strcmp(optimOpt.PenaltyWeight,'user-defined')
         varID = find(w(n_units+1:end));
         for i = varID'
            tw = w(i+n_units);
            if logFlag(i+n_units)
               y = y + tw*(10^(x(i))-varNom(i))^2;
               gy(i) = gy(i) + 2*tw*e10*(10^(x(i))-varNom(i))*10^(x(i));
            else
               y = y + tw*(x(i)-varNom(i))^2;
               gy(i) = gy(i) + 2*tw*(x(i)-varNom(i));
            end
         end
      end
   end

   function [y,gy] = min_1N(x)
      y = 0;
      gy = zeros(nVar+nMD,1);
      varID = find(w(n_units+1:end));
      for i = varID'
         tw = w(n_units+i);
         if logFlag(i+n_units)
            y = y + tw*exp(x(i));
            gy(i) = tw*exp(x(i));
         else
            y = y + tw*abs(x(i));
            gy(i) = tw*sign(x(i));
         end
      end
   end

   function [y,gy] = min_1N2(x)
      y = 0;
      gy = zeros(nVar+nMD,1);
      varID = find(w(n_units+1:end));
      for i = varID'
         tw = w(n_units+i);
         if logFlag(i+n_units)
            y = y + tw*abs(10^(x(i))-varNom(i));
            gy(i) = tw*sign(10^(x(i))-varNom(i))*10^(x(i))*e10;
         else
            y = y + tw*abs(x(i)-varNom(i));
            gy(i) = tw*sign(x(i)-varNom(i));
         end
      end
   end

   function [c,ceq,g,geq] = neq_F(x)
      c = zeros(2*n_units,1);
      g = zeros(nVar+nMD,2*n_units);
      for i = 1:n_units
         id1 = idall{i};
%          if isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.RQModel')
%             N = Nall{i};
%             D = Dall{i};
%             tmpN = [1;x(id1)]'*N*[1;x(id1)];
%             tmpD = [1;x(id1)]'*D*[1;x(id1)];
%             c(2*i-1,1) = tmpN/tmpD - bds(i);
%             c(2*i,1) =  bds(i,1) - tmpN/tmpD;
%             grad1 = zeros(nVar,1);
%             grad2 = zeros(nVar,1);
%             grad1(id1) = (tmpD*N(2:end,2:end)*x(id1)-tmpN*D(2:end,2:end)*x(id1))/tmpD^2;
%             grad2(id1) = -(tmpD*N(2:end,2:end)*x(id1)-tmpN*D(2:end,2:end)*x(id1))/tmpD^2;
%             g(:,2*i-1) = grad1;
%             g(:,2*i) = grad2;
         if isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.QModel')
            quadCoef = Qall{i};
            c(2*i-1,1) = [1;x(id1)]'*quadCoef*[1;x(id1)];
            g(id1,2*i-1) = 2*quadCoef(2:end,2:end)*x(id1)+2*quadCoef(2:end,1);
         end
      end
      c(2:2:end) = bds(:,1)-c(1:2:end);
      c(1:2:end) = c(1:2:end)-bds(:,2);
      g(:,2:2:end) = -g(:,1:2:end);
      ceq = [];
      geq = [];
   end

end






