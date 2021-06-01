function [xs,status] = collectCWAMsamples(obj,n,x0,opt)
 % generate CWAM samples from the feasible set of obj. 
 %  obj - B2BDC.Dataset.Dataset object
 %  n - number of samples wanted
 %  x0 - starting point
 %  opt - B2BDC option
 
 %  Created: Feb, 2017     Wenyu Li
if obj.FeasibleFlag
   opt.AddFitError = true;
end
nsig = opt.SampleOption.AMparam.NumOfSigma;
sig0 = opt.SampleOption.AMparam.Sigma;
warmstart = opt.SampleOption.AMparam.WarmStart;
nStep = opt.SampleOption.StepInterval;
if isempty(warmstart)
   flag1 = false;
elseif all(isfield(warmstart),{'covariance','x0','N','mean','acc'})
   flag1 = true;
else
   error('The specified WarmStart property has incorrect structure field(s)');
end
epsilon = 1e-4;
ac = 0;
if obj.ModelDiscrepancyFlag
   nVar = obj.Variables.Length+obj.ModelDiscrepancy.Variables.Length;
   H = [obj.Variables.calBound; obj.ModelDiscrepancy.Variables.calBound];
else
   nVar = obj.Variables.Length;
   H = [obj.Variables.calBound];
end
if isempty(sig0)
   opt2 = generateOpt;
   xx = obj.collectSample(min(1e4,100*nVar),[],opt2);
   sig0 = cov(xx);
   sv0 = sqrt(diag(sig0))';
elseif all(size(sig0) == nVar)
   sv0 = sqrt(diag(sig0))';
else
   error('The input initial covariance matrix has a wrong dimension');
end
sd = 2.4;
if isempty(x0)
   x0 = obj.FeasiblePoint;
   if obj.ModelDiscrepancyFlag
      x0 = [x0; obj.ModelDiscrepancy.FeasiblePoint];
   end
end
if length(x0) ~= nVar
   error('Wrong dimension for the starting point')
end
if size(x0,1) > 1
   x0 = x0';
end
xs = zeros(n,nVar);
f = likehood_gauY(obj);
n_init = max(5e2,10*nVar);
xt0 = zeros(n_init, nVar);
tmpX = eye(nVar);
if ~flag1
   tcand = rand(n_init*nVar,1);
   xcand = randn(n_init,nVar).*sv0;
   f0 = f(x0);
   count = 1;
   for i = 1:n_init
      for j = 1:nVar
         f1 = f(x0+xcand(i,j)*tmpX(j,:));
         if exp(f1-f0) > tcand(count)
            x0(j) = x0(j)+xcand(i,j);  % accept
            f0 = f1;
         end
         count = count+1;
      end
      xt0(i,:) = x0;
   end
   [ct,t,mx] = calCov(xt0);
else
   ct = warmstart.covariance;
   x0 = warmstart.x0;
   t = warmstart.N;
   mx = warmstart.mean;
   t0 = t;
   acc0 = warmstart.acc;
end
tcand = rand(n*nVar*nStep,1);
f0 = f(x0);
status.NumberOfSample = n;
status.SampleStep = nStep;
count = 1;
if opt.Display
   h = waitbar(0,'Collecting Component-wise Adaptive Metropolis samples...');
   tic;
   for i = 1:n
      for k = 1:nStep
         xcand = randn(1,nVar).*(sqrt(sd*(ct+epsilon)));
         flag = false;
         for j = 1:nVar
            f1 = f(x0+xcand(j)*tmpX(j,:));
            if exp(f1-f0) > tcand(count)
               x0(j) = x0(j)+xcand(j);     % accept
               f0 = f1;
               flag = true;
            end
            count = count+1;
         end
         if flag
            ac = ac+1;
         end
         [ct,t,mx] = updateCov(ct,t,mx,x0);
      end
      xs(i,:) = x0;
      if mod(i,round(0.05*n)) == 0
         waitbar(i/n,h);
      end
   end
   status.CPUTime = toc;
   delete(h);
else
   tic;
   for i = 1:n
      for k = 1:nStep
         xcand = randn(1,nVar).*(sqrt(sd*(ct+epsilon)));
         flag = false;
         for j = 1:nVar
            f1 = f(x0+xcand(j)*tmpX(j,:));
            if exp(f1-f0) > tcand(count)
               x0(j) = x0(j)+xcand(j);     % accept
               f0 = f1;
               flag = true;
            end
            count = count+1;
         end
         if flag
            ac = ac+1;
         end
         [ct,t,mx] = updateCov(ct,t,mx,x0);
      end
      xs(i,:) = x0;
   end
   status.CPUTime = toc;
end
warmstart.covariance = ct;
warmstart.x0 = x0;
warmstart.N = t;
warmstart.mean = mx;
if ~flag1
   warmstart.acc = ac/n*nStep;
else
   warmstart.acc = (acc0*t0+ac)/(t0+n*nStep);
end
status.WarmStart = warmstart;

   function [c,t,mx] = calCov(x)
      t = size(x,1);
      mx = mean(x);
      c = var(x);
%       c = sd*(c+epsilon);
   end

   function [c,t,mxnew] = updateCov(c,t,mx,x)
      mxnew = (t*mx + x)/(t+1);
      c = c*(t-1)/t+mx.^2+x.^2/t-(t+1)/t*mxnew.^2;
      t = t+1;
   end

   function f = likehood_gauY(ds)
      [idx,Qx] = ds.getQ_RQ_expansion;
      nUnits = ds.Length;
      bds = ds.calBound;
      if opt.AddFitError
         for i1 = 1:nUnits
            bds(i1,:) = bds(i1,:) + ds.DatasetUnits.Values(i1).SurrogateModel.ErrorStats.absMax*[-1 1];
         end
      end
      obs = ds.calObserve;
      sigms = (diff(bds,[],2)).^2/2/nsig^2;
      f = @func;
      function y = func(x)
         if any(x>H(:,2)) || any(x<H(:,1))
            y = -inf;
            return
         end
         y = 0;
         for i2 = 1:nUnits
            y = y-([1 x(idx{i2})]*Qx{i2}*[1 x(idx{i2})]'-obs(i2))^2/sigms(i2);
         end
      end
   end

end