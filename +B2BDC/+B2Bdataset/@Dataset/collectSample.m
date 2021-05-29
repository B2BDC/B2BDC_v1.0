function [xSample, status] = collectSample(obj,N,x0,opt)
   % [xSample, status] = collectSample(obj, N) computes uniformly distributed samples
   % from the (joint) feasible set. The samples and associated sampler information are
   % saved in xSample and status, respectively. In this case, the chain starts at the
   % SCM point and the option variable is the default B2BDC option.

   % [xSample, status] = collectSample(obj, N, x0) computes uniformly distributed samples
   % from the (joint) feasible set. The samples and associated sampler information are
   % saved in xSample and status, respectively. In this case, the option variable is the
   % default B2BDC option.
   
   % [xSample, status] = collectSample(obj, N, x0, opt) computes uniformly distributed
   % samples from the (joint) feasible set. The samples and associated sampler information
   % are saved in xSample and status, respectively. In this case, the option variable is
   % specified by opt.

if nargin < 3
   x0 = [];
end
if nargin < 4
   opt = generateOpt('ConsistencyMeasure','absolute');
end
sMD = opt.SampleOption.SampleMethod;
switch sMD
   case 'Gibbs'
      [xSample,status] = collectHRsamples_CW(obj,N,x0,opt);
   case 'HR'
      [xSample,status] = collectHRsamples(obj,N,x0,opt);
   case 'AM'
      [xSample,status] = collectAMsamples(obj,N,x0,opt);
   case 'CWAM'
      [xSample,status] = collectCWAMsamples(obj,N,x0,opt);
end