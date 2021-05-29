function [xSample,status] = collectHRsamples_CW(obj,N,x0,opt)
%  Created: June 28, 2018     Wenyu Li
%  Modified: Oct 15, 2018     Wenyu Li

flag = 0;
CC1 = 0;
CC2 = 0;
% tolX = 1e-2;
% tolY = 0;
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
% UB = UB - 0.5*tolY*(UB-LB);
% LB = LB + 0.5*tolY*(UB-LB);
idQ = [];
AQ = [];
bQ = [];
nV = nVar+nMD;
for i = 1:nUnit
   tCoef = Qx{i}(2:end,2:end);
   if ~any(tCoef(:))
      tA = zeros(1,nV);
      tA(idx{i}) = 2*Qx{i}(1,2:end);
      AQ = [AQ; tA];
      bQ = [bQ; Qx{i}(1,1)];
      idQ = [idQ;i];
   end
end
Ax = [Ax;AQ;-AQ];
bx = [bx;UB(idQ)-bQ;-LB(idQ)+bQ];
idx(idQ) = [];
Qx(idQ) = [];
UB(idQ) = [];
LB(idQ) = [];
nUnit = length(UB);
nStep = opt.SampleOption.StepInterval;
ts = rand(nStep*N*nV,1);
aFlag = cell(nV,1); % which QOI contains the variable
tmpID = 1:nUnit;
for i1 = 1:nV
   tmpFlag = true(nUnit,1);
   for i2 = 1:nUnit
      iid = find(idx{i2}==i1,1);
      if isempty(iid)
         tmpFlag(i2) = false;
      end
   end
   aFlag{i1} = tmpID(tmpFlag);
end
xx = xInit;
xHR = zeros(N,nV);
nS = 1;
status.NumberOfSample = N;
status.SampleStep = nStep;
if opt.Display
   h = waitbar(0,'Collecting Gibbs samples...');
   tic
   for i = 1:nStep*N
      for j = 1:nV
         [tt,tc1,tc2] = calculateIntersection(xx);
         CC1 = CC1+tc1;
         CC2 = CC2+tc2;
         xx(j) = xx(j) + tt;
      end
      if mod(i,nStep) == 0
         xHR(nS,:) = xx';
         nS = nS+1;
      end
      if mod(nS,round(0.05*N)) == 0
         waitbar(nS/N,h);
      end
   end
   delete(h);
   status.CPUTime = toc;
else
   tic;
   for i = 1:nStep*N
      for j = 1:nV
         [tt,tc1,tc2] = calculateIntersection(xx);
         CC1 = CC1 + tc1;
         CC2 = CC2 + tc2;
         xx(j) = xx(j) + tt;
      end
      if mod(i,nStep) == 0
         xHR(nS,:) = xx';
         nS = nS+1;
      end
   end
   status.CPUTime = toc;
end
status.SearcRatio = CC1/CC2;
xSample = xHR;
% xSample.dimension = [nVar nMD];


   function [tt,count1,count2] = calculateIntersection(xx)
      tb = bx-Ax*xx;
      tv = Ax(:,j);
      Pos = tv>0;
      Neg = tv<0;
      tLP = min(tb(Pos)./tv(Pos));
      tLN = max(tb(Neg)./tv(Neg));
      aa = inf(nUnit,1);
      bb = zeros(nUnit,1);
      cUB = ones(nUnit,1);
      cLB = ones(nUnit,1);
      for ii = aFlag{j}
         Coef = Qx{ii};
         iid = find(idx{ii}==j);
         aa(ii) = Coef(iid+1,iid+1);
         bb(ii) = 2*(Coef(iid+1,2:end)*xx(idx{ii})+Coef(1,iid+1));
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
      f = @(t) QuickEval(t,aFlag{j});
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
      tQ(Lflag,:) = repmat([-inf inf],length(Lflag),2);
      tQ = real(tQ);
      aa = aa(:,1);
      bb = bb(:,1);
      tAll = tQ(:);
      tQP = sort(tAll(tAll>0),'ascend');
      tQP(tQP>=tLP) = [];
      nP = length(tQP);
      tQP(nP+1,1) = tLP;
      if nP>0
         tQPm = 0.5*(tQP(1:end-1)+tQP(2:end));
      else
         tQPm = 0.5*tLP;
      end
      tQN = sort(tAll(tAll<=0),'descend');
      tQN(tQN<=tLN) = [];
      nN = length(tQN);
      tQN(nN+1,1) = tLN;
      if nN>0
         tQNm = 0.5*(tQN(1:end-1)+tQN(2:end));
      else
         tQNm = 0.5*tLN;
      end
      count2 = nP+nN;
      count1 = 0;
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
               count1 = count1+1;
            end
            count1 = count1+1;
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
               count1 = count1+1;
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
               count1 = count1+1;
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
               count1 = count1+1;
            end
            count1 = count1+1;
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
               count1 = count1+1;
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
               count1 = count1+1;
               break
            else
               break
            end
         end
      end    
      tn(1,2) = tp(1,2);
      tcand = [tp(2:end,:); tn];
%       if tolX > 0
%          dt = tolX*diff(tcand,[],2);
%          tcand = tcand+[dt -dt];
%       end
      dt = diff(tcand,[],2);
      tsum = [0;cumsum(dt)/sum(dt)];
      it = find(tsum>ts((i-1)*(nV)+j),1);
      if ~isempty(it)
         frac = (ts((i-1)*(nV)+j)-tsum(it-1))/(tsum(it)-tsum(it-1));
      else
         frac = 0;
         flag = flag+1;
      end
      tt = tcand(it-1,1)+dt(it-1)*frac;

      function y = QuickEval(t,flag)
         nn = length(flag);
         y = repmat(aa(flag)*t^2+bb(flag)*t,2,1)+[cUB(flag); cLB(flag)];
         y(nn+1:end) = -y(nn+1:end);
         y = max(y);
      end
   end
end