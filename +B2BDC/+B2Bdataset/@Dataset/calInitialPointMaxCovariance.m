function xInitial = calInitialPointMaxCovariance(obj,n,opt)

% XINITIAL = CALINITIALPOINTMAXCOVARIANCE(OBJ,N,OPT) calculates n initial
% points in the feasible set that maximize the determinant of sampled 
% covariance matrices

%  Created: Aug 12, 2019

fminopt = optimoptions('fmincon','Display','none','GradObj','off','MaxFunctionEvaluations',1e6,...
   'GradConstr','on','Algorithm','interior-point','MaxIter',1e4,'TolFun',1e-10,'TolCon',1e-10,...
   'StepTolerance',1e-20,'DerivativeCheck','off','FiniteDifferenceType','central');
if nargin < 3
   opt = generateOpt;
end
if n < 2
   error('At least 2 initial points are required')
end
nStart = opt.OptimOption.RandomStart;
nCand = 1e2*n;
nVar = obj.Variables.Length;
H = obj.Variables.calBound;
if obj.ModelDiscrepancyFlag
   nMD = obj.ModelDiscrepancy.Variables.Length;
   H = [H; obj.ModelDiscrepancy.Variables.calBound];
else
   nMD = 0;
end
nV = nVar+nMD;
H = repmat(H,n,1);
xOpt = zeros(nStart,nV*n);
yOpt = zeros(nStart,1);
flags = zeros(nStart,1);
vList = obj.Variables;
if ~isempty(vList.ExtraLinConstraint.A)
   tolerance = 1e-5;
   A0 = vList.ExtraLinConstraint.A;
   Aub = vList.ExtraLinConstraint.UB;
   Alb = vList.ExtraLinConstraint.LB;
   eqTest = (Aub-Alb)./sum(abs(A0),2);
   dA = Aub - Alb;
   ieq = (eqTest <= tolerance);
   if any(ieq)
      a1 = [A0(~ieq,:); -A0(~ieq,:)];
      a1(:,nVar+1:nVar+nMD) = 0;
      b1 = [Aub(~ieq); -Alb(~ieq)];
      aeq = A0(ieq,:);
      aeq(:,nVar+1:nVar+nMD) = 0;
      beq = 0.5*(Alb(ieq)+Aub(ieq));
      A1 = [];
      Aeq = [];
      for i = 1:n
         A1 = blkdiag(A1,a1);
         Aeq = blkdiag(Aeq,aeq);
      end
      B1 = repmat(b1,n,1);
      Beq = repmat(beq,n,1);
   else
      a1 = [A0; -A0];
      a1(:,nVar+1:nVar+nMD) = 0;
      Aeq = [];
      b1 = [Aub; -Alb];
      Beq = [];
      A1 = [];
      for i = 1:n
         A1 = blkdiag(A1,a1);
      end
      B1 = repmat(b1,n,1);
   end
else
   A1 = [];
   B1 = [];
   Aeq = [];
   Beq = [];
end
[idall,Qall,Nall,Dall] = obj.getQ_RQ_expansion;
qBD = obj.calBound;
n_units = obj.Length;
abE = zeros(n_units,1);
if opt.AddFitError
   for i = 1:n_units
      if ~isempty(obj.DatasetUnits.Values(i).SurrogateModel.ErrorStats.absMax)
         abE(i) = obj.DatasetUnits.Values(i).SurrogateModel.ErrorStats.absMax;      
      end
   end
end
qBD = qBD+[-abE,abE];
xx = obj.collectHRsamples_CW(nCand,[],opt);
for i = 1:nStart
   x0 = datasample(xx.x,n,'Replace',false);
   x0 = x0';
   [xOpt(i,:),yOpt(i),flags(i)] = fmincon(@funmin2,x0(:),...
      A1,B1,Aeq,Beq,H(:,1),H(:,2),@neq,fminopt);
end
yOpt(flags<=0) = [];
xOpt(flags<=0,:) = [];
if ~isempty(yOpt)
   [~,id] = min(yOpt);
   xInitial = reshape(xOpt(id,:),nV,n);
   xInitial = xInitial';
else
   xInitial = [];
end



   function [y,gy] = funmin(x)
      xs = reshape(x,nV,n);
      xs = xs';
      mx = mean(xs);
      cx = cov(xs);
      y = -det(cx);
      cf = y*inv(cx);
      gy = zeros(nV*n,1);
      for ii = 1:nV
         for jj = ii:nV
            if jj == ii
               gy(ii:nV:end) = gy(ii:nV:end)+2*(xs(:,ii)-mx(ii))*cf(ii,jj);
            else
               gy(ii:nV:end) = gy(ii:nV:end)+2*(xs(:,jj)-mx(jj))*cf(ii,jj);
               gy(jj:nV:end) = gy(jj:nV:end)+2*(xs(:,ii)-mx(ii))*cf(ii,jj);
            end
         end
      end
   end

   function y = funmin2(x)
      xs = reshape(x,nV,n);
      xs = xs';
      cx = cov(xs);
      yy = sort(eig(cx),'descend');
      y = -mean(log(yy(1:min(n-1,nV))));
   end

   function [c,ceq,g,geq] = neq(x)
      c = zeros(2*n_units*n,1);
      g = zeros(nV*n,2*n_units*n);
      for ii = 1:n
         s1 = (ii-1)*nV;
         s2 = 2*(ii-1)*n_units;
         for jj = 1:n_units
            id1 = idall{jj}+s1;
            if isa(obj.DatasetUnits.Values(jj).SurrogateModel,'B2BDC.B2Bmodels.RQModel')
               N = Nall{jj};
               D = Dall{jj};
               c(2*jj-1+s2) = ([1;x(id1)]'*N*[1;x(id1)])/([1;x(id1)]'*D*[1;x(id1)]);
               g(id1,2*jj-1+s2) = 2*((N(2:end,2:end)*x(id1)+N(2:end,1))*([1;x(id1)]'*D*[1;x(id1)]) - ...
                  (D(2:end,2:end)*x(id1)+D(2:end,1))*([1;x(id1)]'*N*[1;x(id1)]))/...
                  ([1;x(id1)]'*D*[1;x(id1)])/([1;x(id1)]'*D*[1;x(id1)]);
            elseif isa(obj.DatasetUnits.Values(jj).SurrogateModel,'B2BDC.B2Bmodels.QModel')
               quadCoef = Qall{jj};
               c(2*jj-1+s2) = [1;x(id1)]'*quadCoef*[1;x(id1)];
               g(id1,2*jj-1+s2) = 2*quadCoef(2:end,2:end)*x(id1)+2*quadCoef(2:end,1);
            end
         end
         c((2:2:2*n_units)+s2) = qBD(:,1)-c((1:2:2*n_units)+s2);
         c((1:2:2*n_units)+s2) = c((1:2:2*n_units)+s2) - qBD(:,2);
      end
      g(:,2:2:end) = -g(:,1:2:end);
      ceq = [];
      geq = [];
   end

end



