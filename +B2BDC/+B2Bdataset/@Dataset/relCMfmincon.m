function [yin_result,s,xopt,abE,flag] = relCMfmincon(obj,disflag,b2bopt)
% subfunction to calculate CM inner bound with fmincon

%  Created: Oct 7, 2016     Wenyu Li

vars = obj.Variables.Values;
vList = obj.Variables;
n_variable = length(vars);
LB = [vars.LowerBound]';
UB = [vars.UpperBound]';
if ~isempty(vList.ExtraLinConstraint.A)
   tolerance = 1e-5;
   A0 = vList.ExtraLinConstraint.A;
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
   n1 = 0.5*size(A1,1);
   n2 = size(Aeq,1);
   A1 = [zeros(2*n1,1), A1];
   if n2 ~= 0
      Aeq = [zeros(n2,1), Aeq];
   end
else
   A1 = [];
   B1 = [];
   Aeq = [];
   Beq = [];
end
units = obj.DatasetUnits.Values;
n_units = length(units);
abE = zeros(length(units),1);
bds = obj.calBound;
bds = diff(bds,[],2);
if b2bopt.AddFitError
   for j = 1:length(units)
      if ~isempty(units(j).SurrogateModel.ErrorStats.absMax)
         abE(j) = units(j).SurrogateModel.ErrorStats.absMax;      
      end
   end
end
flag = quadratictest(obj);

[idall,Qall,Nall,Dall] = obj.getQ_RQ_expansion;

if obj.ModelDiscrepancyFlag
   MDvar = obj.ModelDiscrepancy.Variables;
   nMD = MDvar.Length;
   MDbd = MDvar.calBound;
   LB = [LB; MDbd(:,1)];
   UB = [UB; MDbd(:,2)];
else
   nMD = 0;
end
xbd = UB-LB;
if ~isempty(A1)
   A1 = [A1 zeros(size(A1,1),nMD)]; 
end
if ~isempty(Aeq)
   Aeq = [Aeq zeros(size(Aeq,1),nMD)];
end

% if flag
%    opt = optimoptions('fmincon','Display','none','GradObj','on',...
%       'GradConstr','on','Algorithm','interior-point','Hessian',...
%       'user-supplied','HessFcn',@hessianfcn,'TolFun',1e-6,'TolCon',1e-6);
% else
   opt = optimoptions('fmincon','Display','none','GradObj','on',...
      'GradConstr','on','Algorithm','interior-point','MaxIter',5000,...
      'TolFun',1e-6,'TolCon',1e-6,'DerivativeCheck','off');
% end
if disflag
   disp('=======================================================');
   disp('Searching inner bound...');
   disp('=======================================================');
end
if isempty(obj.FeasiblePoint)
   nStart = b2bopt.OptimOption.RandomStart;
   xStart = zeros(nStart,n_variable+nMD+1);
   xStart(:,2:n_variable+1) = vList.makeDesignSample(nStart);
   ymin = inf;
   for j = 1:nStart
      [x0,yin_result,exitflag,~,ilam] = fmincon(@funxmin,xStart(j,:)',...
         A1,B1,Aeq,Beq,[-Inf;LB],[1;UB],@neq,opt);
      if exitflag > 0
         ymin = min(ymin,yin_result);
      end
      if ymin < 0
         break
      end
   end
else
   xStart = zeros(1,n_variable+nMD+1);
   xStart(2:n_variable+1) = obj.FeasiblePoint;
   if nMD>0
      xStart(n_variable+1+1:n_variable+1+nMD) = obj.ModelDiscrepancy.FeasiblePoint;
   end
   [x0,ymin,exitflag,~,ilam] = fmincon(@funxmin,xStart',...
      A1,B1,Aeq,Beq,[0;LB],[Inf;UB],@neq,opt);
end
xopt = x0(2:n_variable+1);
if exitflag >0
   yin_result = -ymin;
else
   yin_result = -inf;
end
if yin_result > 0
   if nMD > 0
      obj.ModelDiscrepancy.FeasiblePoint = x0(n_variable+2:n_variable+nMD+1);
   end
end
s.expu = ilam.ineqnonlin(1:2:2*n_units).*bds;
s.expl = ilam.ineqnonlin(2:2:2*n_units).*bds;
s.varu = ilam.upper(2:n_variable+nMD+1).*xbd;
s.varl = ilam.lower(2:n_variable+nMD+1).*xbd;
n3 = 0.5*size(A1,1);
if ~isempty(A1)
   s.linu = zeros(n1+n2,1);
   s.linl = zeros(n1+n2,1);
   s.linu(~ieq) = ilam.ineqlin(1:n3).*dA(~ieq);
   s.linl(~ieq) = ilam.ineqlin(n3+1:2*n3).*dA(~ieq);
else
   s.linu = [];
   s.linl = [];
end



   function [y,gy] = funxmin(x)
      y = -x(1);
      gy = [-1;zeros(n_variable+nMD,1)];
   end

   function [c,ceq,g,geq] = neq(x)
      c = zeros(2*n_units,1);
      g = zeros(n_variable+nMD+1,2*n_units);
      for i = 1:n_units
         l = units(i).LowerBound-abE(i);
         u = units(i).UpperBound+abE(i);
         d = units(i).ObservedValue;
         id = idall{i};
         if isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.RQModel')
            N = Nall{i};
            D = Dall{i};
            c(2*i-1,1) = ([1;x(id+1)]'*N*[1;x(id+1)])/([1;x(id+1)]'*D*[1;x(id+1)])-(1-x(1))*u-x(1)*d;
            c(2*i,1) =  (1-x(1))*l+x(1)*d-([1;x(id+1)]'*N*[1;x(id+1)])/([1;x(id+1)]'*D*[1;x(id+1)]);
            grad1 = zeros(n_variable+1,1);
            grad2 = zeros(n_variable+1,1);
            grad1(1) = u-d;
            grad1(id+1) = 2*((N(2:end,2:end)*x(id+1)+N(2:end,1))*([1;x(id+1)]'*D*[1;x(id+1)]) - ...
               (D(2:end,2:end)*x(id+1)+D(2:end,1))*([1;x(id+1)]'*N*[1;x(id+1)]))/...
               ([1;x(id+1)]'*D*[1;x(id+1)])/([1;x(id+1)]'*D*[1;x(id+1)]);
            grad2(1) = d-l;
            grad2(id+1) = -2*((N(2:end,2:end)*x(id+1)+N(2:end,1))*([1;x(id+1)]'*D*[1;x(id+1)]) - ...
               (D(2:end,2:end)*x(id+1)+D(2:end,1))*([1;x(id+1)]'*N*[1;x(id+1)]))/...
               ([1;x(id+1)]'*D*[1;x(id+1)])/([1;x(id+1)]'*D*[1;x(id+1)]);
            g(:,2*i-1) = grad1;
            g(:,2*i) = grad2;
         elseif isa(units(i).SurrogateModel,'B2BDC.B2Bmodels.QModel')
            quadCoef = Qall{i};
            c(2*i-1,1) = [1;x(id+1)]'*quadCoef*[1;x(id+1)]-(1-x(1))*u-x(1)*d;
            c(2*i,1) =  (1-x(1))*l+x(1)*d-[1;x(id+1)]'*quadCoef*[1;x(id+1)];
%             grad1 = zeros(n_variable+nMD+nPD+1,1);
%             grad2 = zeros(n_variable+nMD+nPD+1,1);
%             grad1(id+1) = 2*quadCoef(2:end,2:end)*x(id+1)+2*quadCoef(2:end,1);
%             grad1(1) = u-d;
%             grad2(id+1) = -2*quadCoef(2:end,2:end)*x(id+1)-2*quadCoef(2:end,1);
%             grad2(1) = d-l;
%             g(:,2*i-1) = grad1;
%             g(:,2*i) = grad2;
            g(id+1,2*i-1) = 2*quadCoef(2:end,2:end)*x(id+1)+2*quadCoef(2:end,1);
            g(id+1,2*i) = -g(id+1,2*i-1);
            g(1,2*i-1) = u-d;
            g(1,2*i) = d-l;
         end
      end
      ceq = [];
      geq = [];
   end
end