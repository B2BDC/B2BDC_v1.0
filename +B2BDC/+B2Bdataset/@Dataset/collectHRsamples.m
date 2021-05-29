function [xSample,status] = collectHRsamples(obj,N,x0,opt)

%  Created: June 28, 2018     Wenyu Li
%  Modified: Oct 15, 2018     Wenyu Li

flag = 0;
if ~quadratictest(obj)
   error('HR sampling is now only available for quadratic models')
end
if obj.FeasibleFlag
   opt.AddFitError = true;
end
if obj.isConsistent(opt)
   if obj.ModelDiscrepancyFlag
      nMD = length(obj.ModelDiscrepancy.FeasiblePoint);
   else
      nMD = 0;
   end
   if isempty(x0)
      xInit = obj.FeasiblePoint;
      if nMD > 0
         xInit = [xInit; obj.ModelDiscrepancy.FeasiblePoint];
      end
   else
      if size(x0,1) == 1
         x0 = x0';
      end
      [~,xInit] = obj.isFeasiblePoint(x0');
      if isempty(xInit)
         error('The input starting point is not feasible');
      else
         xInit = xInit';
      end
   end
else
   errordlg('Infeasible problem: Dataset cannot be shown to be Consistent');
   error('collectSamples:Inconclusive',...
      'Dataset cannot be shown to be Consistent');
end
% tolX = 0;
% tolY = 1e-4;
[idx,Qx] = obj.getQ_RQ_expansion;
bounds = obj.calBound;
LB = bounds(:,1);
UB = bounds(:,2);
sVar = obj.Variables;
bounds = sVar.calBound;
if nMD > 0
   bounds = [bounds; obj.ModelDiscrepancy.Variables.calBound];
end
nVar = sVar.Length;
Ax = [eye(nVar+nMD); -eye(nVar+nMD)];
bx = [bounds(:,2); -bounds(:,1)];
if ~isempty(sVar.ExtraLinConstraint.A)
   AA = [sVar.ExtraLinConstraint.A zeros(size(sVar.ExtraLinConstraint.A,1),nMD)];
   Ax = [Ax; AA; -AA];
   bx = [bx; sVar.ExtraLinConstraint.UB; -sVar.ExtraLinConstraint.LB];
end
units = obj.DatasetUnits.Values;
nUnit = length(units);
if obj.FeasibleFlag
   for i = 1:nUnit
      if ~isempty(units(i).SurrogateModel.ErrorStats)
         LB(i) = LB(i)-units(i).SurrogateModel.ErrorStats.absMax;
         UB(i) = UB(i)+units(i).SurrogateModel.ErrorStats.absMax;
      end
   end
end
idQ = [];
AQ = [];
nV = nVar+nMD;
for i = 1:nUnit
   tCoef = Qx{i}(2:end,2:end);
   if ~any(tCoef(:))
      tA = zeros(1,nV);
      tA(idx{i}) = 2*tCoef(1,2:end);
      AQ = [AQ; tA];
      idQ = [idQ;i];
   end
end
Ax = [Ax;AQ;-AQ];
bx = [bx;UB(idQ);-LB(idQ)];
idx(idQ) = [];
Qx(idQ) = [];
UB(idQ) = [];
LB(idQ) = [];
nUnit = length(UB);
nStep = opt.SampleOption.StepInterval;
mu = zeros(1,nV);
% if nargin < 5
sigma = eye(nV);
d = mvnrnd(mu,sigma,nStep*N);
d = d';
% else
%    sigma = PC.sigma;
%    vPC = PC.vv;
%    d = mvnrnd(mu,sigma,nStep*N);
%    d = vPC*d';
% end
% d = randn(nV,nStep*N);
d = normc(d);
% ts = (1-0.5*tolerance)*rand(nStep*N,1)+0.5*tolerance;
ts = rand(nStep*N,1);
Ad = Ax*d;
xx = xInit;
xHR = zeros(N,nV);
nS = 1;
status.NumberOfSample = N;
status.SampleStep = nStep;
CC1 = 0;
CC2 = 0;
if opt.Display
   tic;
   h = waitbar(0,'Collecting hit and run samples...');
   for i = 1:nStep*N
      [tt,tcc1,tcc2] = calculateIntersection(d(:,i),xx);
      CC1 = CC1+tcc1;
      CC2 = CC2+tcc2;
      xx = xx + tt*d(:,i);
      if mod(i,nStep) == 0
         xHR(nS,:) = xx';
         nS = nS+1;
      end
      if mod(nS,round(0.05*N)) == 0
         waitbar(nS/N,h);
      end
   end
   status.CPUTime = toc;
else
   tic;
   for i = 1:nStep*N
      [tt,tcc1,tcc2] = calculateIntersection(d(:,i),xx);
      CC1 = CC1+tcc1;
      CC2 = CC2+tcc2;
      xx = xx + tt*d(:,i);
      if mod(i,nStep) == 0
         xHR(nS,:) = xx';
         nS = nS+1;
      end
   end
   status.CPUTime = toc;
end
status.SearcRatio = CC1/CC2;
% if sum(obj.isFeasiblePoint(xHR)) < N
%    disp('Errors');
% end
xSample = xHR;
% xSample.dimension = [nVar nMD];

   function [tt,count,count2] = calculateIntersection(d,xx)
      tb = bx-Ax*xx;
      tv = Ad(:,i);
      Pos = tv>0;
      Neg = tv<0;
      tLP = min(tb(Pos)./tv(Pos));
      tLN = max(tb(Neg)./tv(Neg));
      aa = zeros(nUnit,1);
      bb = zeros(nUnit,1);
      cUB = zeros(nUnit,1);
      cLB = zeros(nUnit,1);
      for ii = 1:nUnit
         Coef = Qx{ii};
         vv = d(idx{ii});
         aa(ii) = vv'*Coef(2:end,2:end)*vv;
         bb(ii) = 2*(vv'*Coef(2:end,2:end)*xx(idx{ii})+Coef(1,2:end)*vv);
         cc = [1;xx(idx{ii})]'*Coef*[1;xx(idx{ii})];
         cUB(ii) = cc-UB(ii);
         cLB(ii) = cc-LB(ii);
      end
      Lflag = find(aa==0);
      if ~isempty(Lflag)
         tL = [-cUB(Lflag)./bb(Lflag); -cLB(Lflag)./bb(Lflag)];
         tLP = min(tLP,min(tL(tL>0)));
         tLN = max(tLN,max(tL(tL<0)));
      end
      f = @QuickEval;
      delta_UB = bb.^2-4*aa.*cUB;
      delta_LB = bb.^2-4*aa.*cLB;
      flag_UB = find(delta_UB <= 0);
      flag_LB = find(delta_LB <= 0);
      aa = repmat(aa,1,4);
      bb = repmat(bb,1,4);
      cc = repmat([-1,1],nUnit,2).*sqrt([delta_UB delta_UB delta_LB delta_LB]);
      tQ = 0.5*(-bb + cc)./aa;
      tQ(flag_UB,1:2) = repmat([-inf inf],length(flag_UB),1);
      tQ(flag_LB,3:4) = repmat([-inf inf],length(flag_LB),1);
      tQ = real(tQ);
      aa = aa(:,1);
      bb = bb(:,1);
      tAll = tQ(:);
      tQP = sort(tAll(tAll>0),'ascend');
      tQP(tQP>=tLP) = [];
      nP = length(tQP);
      tQP(end+1,1) = tLP;
      if nP>0
         tQPm = 0.5*(tQP(1:end-1)+tQP(2:end));
      else
         tQPm = 0.5*tLP;
      end
      tQN = sort(tAll(tAll<=0),'descend');
      tQN(tQN<=tLN) = [];
      nN = length(tQN);
      tQN(end+1,1) = tLN;
      if nN>0
         tQNm = 0.5*(tQN(1:end-1)+tQN(2:end));
      else
         tQNm = 0.5*tLN;
      end
      count = 0;
      count2 = nP+nN;
      % positive t
      tp = [0 0];
      idp = 1;
      if nP == 0
         tp(end,2) = tQP;
      else
         while idp <= nP
            % search for infeasible right end
            while idp <= nP && f(tQPm(idp))<0
               idp = idp+1;
               count = count+1;
            end
            count = count+1;
            if idp > nP
               tp(end,2) = tQP(end);
               break
            else
               tp(end,2) = tQP(idp);
            end
            % find next candidate feasible left end point
            [ti,tj] = find(tQ==tp(end,2));
            if aa(ti) > 0
               if tj == 3
                  if tQ(ti,1) > tQ(ti,tj) && tQ(ti,1) < tQ(ti,4) && ~isempty(find(tQP == tQ(ti,1),1))
                     idp = find(tQP == tQ(ti,1),1);
                  elseif ~isempty(find(tQP == tQ(ti,4),1))
                     idp = find(tQP == tQ(ti,4),1);
                  else
                     break;
                  end
               else
                  break
               end
            else
               if tj == 1
                  if tQ(ti,3) > tQ(ti,tj) && tQ(ti,3) < tQ(ti,2) && ~isempty(find(tQP == tQ(ti,3),1))
                     idp = find(tQP == tQ(ti,3),1);
                  elseif ~isempty(find(tQP == tQ(ti,2),1))
                     idp = find(tQP == tQ(ti,2),1);
                  else
                     break;
                  end
               else
                  break;
               end
            end
            % find a feasible left end point
            id0 = idp;
            for jj = idp:nP-1
               yy = f(tQPm(jj));
               count = count+1;
               if yy<0
                  tp(end+1,1) = tQP(jj);
                  idp = jj+1;
                  break
               end
            end
            if idp > id0
               continue
            elseif f(tQPm(end))<0
               tp(end+1,:) = tQP(end-1:end);
               count = count+1;
               break
            else
               break
            end
         end
      end
      
      % negative t
      tn = [0 0];
      idn = 1;
      if nN == 0
         tn(end,1) = tQN;
      else
         while idn <= nN
            % search for infeasible left end
            while idn <= nN && f(tQNm(idn))<0
               idn = idn+1;
               count = count+1;
            end
            count = count+1;
            if idn > nN
               tn(end,1) = tQN(end);
               break
            else
               tn(end,1) = tQN(idn);
            end
            % find next candidate feasible right end point
            [ti,tj] = find(tQ==tn(end,1));
            if aa(ti) > 0
               if tj == 4
                  if tQ(ti,2) < tQ(ti,tj) && tQ(ti,2) > tQ(ti,3) && ~isempty(find(tQN == tQ(ti,2),1))
                     idn = find(tQN == tQ(ti,2),1);
                  elseif ~isempty(find(tQN == tQ(ti,3),1))
                     idn = find(tQN == tQ(ti,3),1);
                  else
                     break;
                  end
               else
                  break;
               end
            else
               if tj == 2
                  if tQ(ti,4) < tQ(ti,tj) && tQ(ti,4) > tQ(ti,1) && ~isempty(find(tQN == tQ(ti,4),1))
                     idn = find(tQN == tQ(ti,4),1);
                  elseif ~isempty(find(tQP == tQ(ti,1),1))
                     idn = find(tQP == tQ(ti,1),1);
                  else
                     break
                  end
               else
                  break
               end
            end
            % find a feasible right end point
            id0 = idn;
            for jj = idn:nN-1
               yy = f(tQNm(jj));
               count = count+1;
               if yy<0
                  tn(end+1,2) = tQN(jj);
                  idn = jj+1;
                  break
               end
            end
            if idn > id0
               continue
            elseif f(tQNm(end))<0
               tn(end+1,:) = flipud(tQN(end-1:end));
               count = count+1;
               break
            else
               break
            end
         end
      end
      tn(1,2) = tp(1,2);
      tcand = [tp(2:end,:); tn];
      dt = diff(tcand,[],2);
      tsum = [0;cumsum(dt)/sum(dt)];
      it = find(tsum>ts(i),1);
      if ~isempty(it)
         frac = (ts(i)-tsum(it-1))/(tsum(it)-tsum(it-1));
         % not too close to the boundary
%          if frac < 0.5*tolX
%             frac = 0.5*tolX;
%          elseif frac > 1-0.5*tolX
%             frac = 1-0.5*tolX;
%          end
      else
         frac = 0;
         flag = flag+1;
      end
      tt = tcand(it-1,1)+dt(it-1)*frac;
        
         function y = QuickEval(t)
            y = repmat(aa*t^2+bb*t,2,1) + [cUB; cLB];
            y(nUnit+1:end) = -y(nUnit+1:end);
            y = max(y);
         end
   end
end