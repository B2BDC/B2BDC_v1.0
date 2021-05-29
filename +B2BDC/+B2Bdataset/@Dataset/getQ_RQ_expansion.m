function [idall,Qall,Nall,Dall] = getQ_RQ_expansion(obj,Qpred)
% [IDALL,QALL,NALL,DALL] = GETQ_RQ_EXPANSION(OBJ) calculates the corresponding
% variable index, quadratic coefficient matrix, rational quadratic numerator 
% coefficient matrix and rational quadratic denumenator coefficient matrix for
% consistency measure or prediction calculation.

% Created: Sep 9, 2018     Wenyu Li


allVarnames = obj.VarNames;
n_variable = length(allVarnames);
units = obj.DatasetUnits.Values;
n_units = length(units);
if nargin > 1
   n_units = n_units+1;
   model0 = Qpred.Model;
   predGroup = Qpred.Correction.GroupIndex;
   sv0 = Qpred.Correction.Value;
end
idall = cell(n_units,1);
Nall = cell(n_units,1);
Dall = cell(n_units,1);
Qall = cell(n_units,1);
if obj.ModelDiscrepancyFlag
   MDinfo = obj.ModelDiscrepancy;
   MDBasis = MDinfo.Basis;
   MDnames = {MDinfo.Variables.Values.Name}';
   nmd = MDinfo.CorrectionDimension;
   if nargin > 1 && predGroup > 0
      fMD = MDinfo.BasisFunction;
      MDBasis{end+1} = fMD{predGroup}(sv0);
   end
   GroupIndex = MDinfo.GroupIndex;
   for j = 1:n_units
      if j==n_units && nargin > 1
         tmodel = model0;
         [~,~,id1] = intersect(tmodel.VarNames,[allVarnames; MDnames],'stable');
         id3 = id1(id1>n_variable)-n_variable;
         id1 = id1(id1<=n_variable);
      else
         tmodel = units(j).SurrogateModel;
         [~,~,id1] = intersect(tmodel.VarNames,allVarnames,'stable');
      end
      if j==n_units && nargin>1
         gid = predGroup;
      else
         gid = GroupIndex(j);
      end
      if gid > 0
         id2 = (sum(nmd(1:gid-1))+1:sum(nmd(1:gid)))';
      else
         id2 = [];
      end
      
      if isa(tmodel,'B2BDC.B2Bmodels.RQModel')
         error('Model discrepancy correction is currently not available for rational quadratic models');
      elseif isa(tmodel,'B2BDC.B2Bmodels.QModel')
         if j == n_units && nargin > 1 && ~isempty(id3)
            if isempty(id2)
               Coef = tmodel.CoefMatrix;
               idall{j} = [id1; id3+n_variable];
            else
               id4 = sort(unique([id2;id3]));
               idall{j} = [id1;id4+n_variable];
               nc1 = length(id4);
               tCoef = tmodel.CoefMatrix;
               nm = size(tCoef,1);
               Coef = zeros(nm+nc1+1);
               [~,~,id3] = intersect(id3,id4,'stable');
               [~,~,id2] = intersect(id2,id4,'stable');
               Coef([0:nm, id3'+nm]+1,[0:nm, id3'+nm]+1) = tCoef;
               Coef(1,id2+nm+1) = Coef(1,id2+1)+0.5*MDBasis{j};
               Coef(id2+nm+1,1) = Coef(id2+1,1)+0.5*MDBasis{j};
            end
         else
            nc1 = length(id2);
            tCoef = tmodel.CoefMatrix;
            nm = size(tCoef,1);
            Coef = zeros(nm+nc1);
            Coef(1:nm,1:nm) = tCoef;
            if nc1 > 0
               Coef(1,nm+1:nm+nc1) = 0.5*MDBasis{j};
               Coef(nm+1:nm+nc1,1) = 0.5*MDBasis{j};
            end
            idall{j} = [id1;id2+n_variable];
         end
         Qall{j} = Coef;
      end
   end
else
   for j = 1:n_units
      if j == n_units && nargin > 1
         tmodel = model0;
      else
         tmodel = units(j).SurrogateModel;
      end
      [~,~,idall{j}] = intersect(tmodel.VarNames,allVarnames,'stable');
      if isa(tmodel,'B2BDC.B2Bmodels.RQModel')
         Nall{j} = tmodel.Numerator;
         Dall{j} = tmodel.Denominator;
      elseif isa(tmodel,'B2BDC.B2Bmodels.QModel')
         Qall{j} = tmodel.CoefMatrix;
      end
   end
end