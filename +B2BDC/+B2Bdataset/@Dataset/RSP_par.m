function [xVal,eff] = RSP_par(obj,N,xS,Dinfo,opt)
% sub function for feasible set uniform sampling

%  Created: Oct 18, 2016     Wenyu Li


p = parpool('IdleTimeout', 120);
npool = p.NumWorkers;
xVal = [];
if npool >= 4
   nt1 = npool;
else
   nt1 = 2*npool;
end
th = opt.SampleOption.RejectionTol;
da = opt.SampleOption.DirectionAdaption;
if da > 0
   daFlag = true;
else
   daFlag = false;
end
sopt = opt.SampleOption;
s1 = sopt.UncertaintyEstimation;
vars = obj.Variables;
nVar = vars.Length;
nCut = sopt.ExtraCut;
nPC = nVar-opt.SampleOption.TruncatedPC;
if ~isempty(sopt.DataStorePath) && isdir(sopt.DataStorePath)
   savePath = sopt.DataStorePath;
   fid = fopen(fullfile(savePath,'Computation performance.txt'),'w');
else
   savePath = [];
   fid = -1;
end
if nPC < 1
   nPC = 1;
end
if nPC < nVar ||...
      ((strcmp(s1,'Sample') || strcmp(s1,'Truncation')) && nCut >0)
   if isempty(xS)
      nS = min(10^6,10*nVar^3);
      xS = obj.collectSamples(nS,[],opt);
   else
      nS = size(xS,1);
   end
end
if nPC < nVar
   xAve = mean(xS);
   if ~obj.isFeasiblePoint(xAve)
      Xdd = xS - repmat(xAve,nS,1);
      xdd = sqrt(sum(Xdd.^2,2));
      [~,imin] = min(xdd);
      xAve = xS(imin,:);
   end
   xC = xS - repmat(xAve,nS,1);
   [V,DD] = eig((xC' * xC)/nS);
   Ddiag = diag(DD);
   [~,id] = sort(Ddiag,'descend');
   vv = V(:,id);
   y0 = xAve*vv;
   sampledPC.V = vv;
   sampledPC.D = Ddiag(id);
   sampledPC.x0 = xAve;
   if ~isempty(savePath)
      save(fullfile(savePath,'PCinfo'),'sampledPC');
   end
else
   vv = eye(nVar);
end
if ~isempty(Dinfo)
   if size(Dinfo.D,1) == nPC
      A = Dinfo.D;
      uq = Dinfo.uq;
   else
      disp('The input D matrix should have nPC rows');
      A = [];
      uq = [];
   end
else
   A = [];
   uq = [];
end
if nPC < nVar
   projectDS = obj.projectDSonActiveSubspace(sampledPC,nPC);
   vars = projectDS.Variables;
   uqv = zeros(nPC,2);
   switch sopt.UncertaintyEstimation
      case 'Outer'
         opt.Prediction = 'outer';
         for i = 1:nPC
            tv = zeros(nPC+1,1);
            tv(i+1) = 1;
            tM = generateModel(tv,vars);
            preQ = projectDS.predictQOI(tM,opt);
            uqv(i,:) = [preQ.min preQ.max];
         end
      case 'Inner'
         opt.Prediction = 'inner';
         for i = 1:nPC
            tv = zeros(nPC+1,1);
            tv(i+1) = 1;
            tM = generateModel(tv,vars);
            preQ = projectDS.predictQOI(tM,opt);
            uqv(i,:) = [preQ.min preQ.max];
         end
      case 'Sample'
         xx = xS-repmat(mean(xS),size(xS,1),1);
         xC = xx*vv(:,1:nPC);
         uqv = [min(xC)' max(xC)'];
      case 'Truncation'
         sf = 1-sopt.UncertaintyTruncation;
         opt.Prediction = 'inner';
         for i = 1:nCut
            tv = zeros(nPC+1,1);
            tv(i+1) = 1;
            tM = generateModel(tv,vars);
            preQ = projectDS.predictQOI(tM,opt);
            fmin = preQ.min;
            fmax = preQ.max;
            uqv(i,:) = sf*[fmin, fmax];
         end
   end
   sVar = generateVar([],uqv);
else
   sVar = obj.Variables;
end
if nCut > 0
   tic
   if nPC < nVar
      subDS = projectDS;
   else
      subDS = obj;
   end
   uq2 = zeros(nCut,2);
   D = randn(nPC,nCut);
   D = normc(D);
   switch sopt.UncertaintyEstimation
      case 'Outer'
         opt.Prediction = 'outer';
         parfor i = 1:nCut
            tv = [0; D(:,i)];
            tM = generateModel(tv,sVar);
            preQ = subDS.predictQOI(tM,opt);
            uq2(i,:) = [preQ.min preQ.max];
         end
      case 'Inner'
         opt.Prediction = 'inner';
         parfor i = 1:nCut
            tv = [0; D(:,i)];
            tM = generateModel(tv,sVar);
            preQ = subDS.predictQOI(tM,opt);
            uq2(i,:) = [preQ.min preQ.max];
         end
      case 'Sample'
         parfor i = 1:nCut
            f = xC*D(:,i);
            uq2(i,:) = [min(f) max(f)];
         end
      case 'Truncation'
         opt.Prediction = 'inner';
         parfor i = 1:nCut
            tv = [0; D(:,i)];
            tM = generateModel(tv,sVar);
            preQ = subDS.predictQOI(tM,opt);
            fmin = preQ.min;
            fmax = preQ.max;
            uq2(i,:) = sf*[fmin, fmax];
         end
   end
   tt = toc;
   if fid > 0
      fprintf(fid,'The CPU time to calcualte the polytope with %d directions is: %5.3E \n',nCut,tt);
   end
else
   D = [];
   uq2 = [];
end
A = [A D]';
uq = [uq;uq2];
if size(A,1) > 0
   sVar = sVar.addLinearConstraint(A,uq(:,1),uq(:,2));
   if ~isempty(savePath)
      save(fullfile(savePath,'polytopeInfo'),'sVar');
   end
end
ns = sopt.BatchMaxSample;
tic;
xx = [];
ss = generateSampleOpt('StepInterval',100);
while any(daFlag)
   nt1 = ceil(da*N/th/ns);
   xApp = [];
   xx = sVar.collectSamples(nt1,xx(end,:),ss);
   xCand = cell(nt1,1);
   parfor i = 1:nt1
      xcand = sVar.collectSamples(ns,xx(i,:),opt.SampleOption);
      if nPC < nVar
         xtmp = repmat(xAve,ns,1)+xcand*vv(:,1:nPC)';
      else
         xtmp = xcand;
      end
      iF = obj.isFeasiblePoint(xtmp);
      xCand{i} = xtmp(iF,:);
   end
   for i = 1:nt1
      xApp = [xApp; xCand{i}];
   end
   if size(xApp,1)<0.9*da*N
      eff = size(xApp,1)/nt1/ns;
      disp(['Numerical efficiency of current sampling method is ' num2str(eff)]);
      return;
   end
   [sVar,daFlag] = sVar.DirectionalAdapt(xApp,daFlag);
end
tt = toc;
if fid > 0
   fprintf(fid,'The CPU time to update polytope directional facets is: %5.3E \n',tt);
end
tic;
nt1 = 2*npool;
xCand = cell(nt1,1);
xx = sVar.collectSamples(nt1,xx(end,:),ss);
parfor i = 1:nt1
   xcand = sVar.collectSamples(ns,xx(i,:),opt.SampleOption);
   if nPC < nVar
      xtmp = repmat(xAve,ns,1)+xcand*vv(:,1:nPC)';
   else
      xtmp = xcand;
   end
   iF = obj.isFeasiblePoint(xtmp);
   xCand{i} = xtmp(iF,:);
end
for i = 1:nt1
   xVal = [xVal;xCand{i}];
   if ~isempty(savePath)
      save(fullfile(savePath,'xSample'),'xVal');
   end
end
tt = toc;
n1 = size(xVal,1);
eff = n1/nt1/ns;
if eff < th
   disp(['Numerical efficiency of current sampling method is ' num2str(eff)])
   delete(p)
   return;
elseif n1 >= N
   xVal = xVal(randperm(n1,N),:);
   if ~isempty(savePath)
      save(fullfile(savePath,'efficiency'),'eff');
   end
   delete(p)
   if fid > 0
      fclose(fid);
   end
   return;
else
   if ~isempty(savePath)
      save(fullfile(savePath,'efficiency'),'eff');
   end
   nEst = ceil(1.2*(N-n1)/eff/ns);
   nr = 0.5*ceil(nEst/npool);
   if fid > 0
      fprintf('Current progress: %4.2f%%, estimated time left: %5.3E \n',n1/N,tt*nr);
   end
   xCand = cell(nEst,1);
   xx = sVar.collectSamples(nEst,xx(end,:),ss);
   parfor i = 1:nEst
      xcand = sVar.collectSamples(ns,xx(i,:),opt.SampleOption);
      if nPC < nVar
         xtmp = repmat(xAve,ns,1)+xcand*vv(:,1:nPC)';
      else
         xtmp = xcand;
      end
      iF = obj.isFeasiblePoint(xtmp);
      xCand{i} = xtmp(iF,:);
   end
   for i = 1:nEst
      xVal = [xVal; xCand{i}];
   end
   n1 = size(xVal,1);
   eff = n1/nt1/ns;
   if ~isempty(savePath)
      save(fullfile(savePath,'efficiency'),'eff');
   end
   if ~isempty(savePath)
      save(fullfile(savePath,'xSample'),'xVal');
   end
   if n1>N
      xVal = xVal(randperm(n1,N),:);
   end
   if fid > 0
      fclose(fid);
   end
   delete(p);
end

