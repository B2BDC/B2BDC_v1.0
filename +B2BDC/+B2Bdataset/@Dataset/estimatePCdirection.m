function [xFea,v,sv,cumsv,xE] = estimatePCdirection(obj,opt,nS)

% estimate PC directions based on uniform direction and directional line
% segment distance.

% xFea: the center point (consistency measure point)
% v: estimated principal direction, stored in column vectors in v
% sv: estimated length scale associated with each principal direction
% cumsv: cumulative squared length scale, representing TEV
% xE: the extreme points along each principal direction
% nS (optional): number of candidate directions in each step

%   Created: Jan 11, 2019     Wenyu Li


if nargin < 2
   opt = generateOpt('Display',false,'Prediction','inner','LocalStart',5);
end
if ~obj.isConsistent(opt)
   error('The dataset is inconsistent');
else
   xFea = obj.FeasiblePoint;
end
nVar = obj.Variables.Length;
if nargin < 3
   nS = 5*nVar;
end
if obj.ModelDiscrepancyFlag
   nMD = obj.ModelDiscrepancy.Variables.Length;
   xFea = [xFea; obj.ModelDiscrepancy.FeasiblePoint];
else
   nMD = 0;
end
if obj.ParameterDiscrepancyFlag
   nPD = obj.ParameterDiscrepancy.Variables.Length;
   xFea = [xFea; obj.ParameterDiscrepancy.FeasiblePoint];
else
   nPD = 0;
end
v = zeros(nVar+nMD+nPD);
sv = zeros(nVar+nMD+nPD,1);
xE = zeros(2*(nVar+nMD+nPD),nVar+nMD+nPD);
[idx,Qx,~,~,APD,bPD] = obj.getQ_RQ_expansion;
bounds = obj.calBound;
LB = bounds(:,1);
UB = bounds(:,2);
sVar = obj.Variables;
bounds = sVar.calBound;
nVar = sVar.Length;
Ax = [eye(nVar); -eye(nVar)];
Ax = [Ax zeros(2*nVar,nMD+nPD)];
bx = [bounds(:,2); -bounds(:,1)];
if ~isempty(sVar.ExtraLinConstraint.A)
   AA = [sVar.ExtraLinConstraint.A zeros(size(sVar.ExtraLinConstraint.A,1),nMD+nPD)];
   Ax = [Ax; AA; -AA];
   bx = [bx; sVar.ExtraLinConstraint.UB; -sVar.ExtraLinConstraint.LB];
end
Ax = [Ax; APD];
bx = [bx; bPD];
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
nD = randn(nVar+nMD+nPD,nS);
nD = normc(nD);
Ad = Ax*nD;
LL = zeros(nS,1);
tmins = zeros(nS,1);
tmaxs = zeros(nS,1);
for i = 1:nS
   [LL(i),tmins(i),tmaxs(i)] = calculateIntersection(nD(:,i),xFea);
end
[sv(1),id] = max(LL);
v(:,1) = nD(:,id);
xE(1,:) = xFea+nD(:,id)*tmins(id);
xE(nVar+nMD+nPD+1,:) = xFea+nD(:,id)*tmaxs(id);
for j = 2:nVar
   Z = null(v(:,1:j-1)');
   nD = randn(nVar+nMD+nPD-j+1,nS);
   nD = normc(Z*nD);
   Ad = Ax*nD;
   LL = zeros(nS,1);
   tmins = zeros(nS,1);
   tmaxs = zeros(nS,1);
   for i = 1:nS
      [LL(i),tmins(i),tmaxs(i)] = calculateIntersection(nD(:,i),xFea);
   end
   [sv(j),id] = max(LL);
   v(:,j) = nD(:,id);
   xE(j,:) = xFea+nD(:,id)*tmins(id);
   xE(nVar+nMD+nPD+j,:) = xFea+nD(:,id)*tmaxs(id);
end
[sv,id] = sort(sv,'descend');
xE(id,:) = xE(1:nVar,:);
xE(id+nVar,:) = xE(nVar+1:end,:);
v = v(:,id);
cumsv = cumsum(sv.^2);
cumsv = cumsv/cumsv(end);

   function [tt,tmin,tmax] = calculateIntersection(d,xx)
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
      % positive t
      tp = [0 0];
      idp = 1;
      if nP == 0
         tp(end,2) = tQP;
      else
         while idp <= nP
            % search for infeasible right end
            while idp <= nP && all(f(tQPm(idp))<0)
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
               if all(yy<0)
                  tp(end+1,1) = tQP(jj);
                  idp = jj+1;
                  break
               end
            end
            if idp > id0
               continue
            elseif all(f(tQPm(end))<0)
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
            while idn <= nN && all(f(tQNm(idn))<0)
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
               if all(yy<0)
                  tn(end+1,2) = tQN(jj);
                  idn = jj+1;
                  break
               end
            end
            if idn > id0
               continue
            elseif all(f(tQNm(end))<0)
               tn(end+1,:) = flipud(tQN(end-1:end));
               count = count+1;
               break
            else
               break
            end
         end
      end
      tcand = [tp; tn];
      dt = diff(tcand,[],2);
      tt = sum(dt);
      tmin = tn(end,1) + 0.01*diff(tn(end,:));
      tmax = tp(end,2) - 0.01*diff(tp(end,:));
      
      function y = QuickEval(t)
         y = zeros(2*nUnit,1);
         y(1:nUnit) = aa*t^2+bb*t+cUB;
         y(nUnit+1:end) = -(aa*t^2+bb*t+cLB);
      end
   end
end