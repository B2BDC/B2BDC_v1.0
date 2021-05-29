function dsUnitList = checkSelfConsistency(obj,opt)
%  dsUnitList = checkSelfConsistency(obj) calculates self-inconsistent dataset units in the
%  dataset. The option variable is the default B2BDC option object.
%
%  dsUnitList = checkSelfConsistency(obj,opt) calculates self-inconsistent dataset units
%  in the dataset. The option variable is specified by the input opt

% Created: Oct 20, 2015      Wenyu Li

if nargin < 2
   opt = generateOpt;
   opt.Display = false;
   opt.ExtraLinFraction = -1;
end
dsUnitList = [];
n = obj.Length;
if opt.Display
   h = waitbar(0,'Start to check dataset unit self-consistency');
end
for i = 1:n
   dsTest = B2BDC.B2Bdataset.Dataset;
   dsUnit = obj.DatasetUnits.Values(i);
   dsTest.addDSunit(dsUnit);
   dsTest.isConsistent(opt);
   if dsTest.ConsistencyMeasure(2) <= 0
      dsUnitList = [dsUnitList; dsUnit];
   end
   if opt.Display
      waitbar(i/n,h,['Checking dataset unit self-consistency ' num2str(i) ' / ' num2str(n)]);
   end
end
%   if opt.Display
%      waitbar(1,h,'Deleting self-inconsistent dataset units from dataset');
%   end
if opt.Display
   close(h);
end