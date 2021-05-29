function y = eval_with_discrepancy(obj,X,DSIdx)
%  y = eval with discrepancy(obj, X) evaluates all QOI outputs at input design points
%  specified by X. Each row of X corresponds to a design point including both model
%  parameter and discrepancy-function coefficients. Each column of y corresponds to a
%  QOI.
%
%  y = eval with discrepancy(obj, X, DSIdx) evaluates some QOI outputs, specified by
%  DSIdx, at input design points specified by X. Each row of X corresponds to a design
%  point including both model parameter and discrepancy-function coefficients. Each
%  column of y corresponds to a QOI.

%  Created: Sep 12, 2018     Wenyu Li
[idall,Qall] = obj.getQ_RQ_expansion;
if nargin > 2
   if iscell(DSIdx)
      allDS = {obj.DatasetUnits.Values.Name}';
      [~,~,idQOI] = intersect(DSIdx,allDS,'stable');
   else
      [~,~,idQOI] = intersect(DSIdx,1:obj.Length,'stable');
   end
   if length(idQOI) ~= length(DSIdx)
      error('Wrong specified QOI index');
   end
else
   idQOI = 1:obj.Length;
end
nSample = size(X,1);
nQOI = length(idQOI);
y = zeros(nSample,nQOI);
for i = 1:nQOI
   Q0 = Qall{idQOI(i)};
   id0 = idall{idQOI(i)};
   xx = B2BDC.Fitting.expandBasis(X(:,id0));
   vec = B2BDC.Fitting.coef2vec(Q0);
   y(:,i) = xx*vec;
end