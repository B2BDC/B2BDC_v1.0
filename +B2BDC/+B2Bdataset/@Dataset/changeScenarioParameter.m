function changeScenarioParameter(obj,Svalue,Sname,QOIindex)
%  changeScenarioParameter(obj,Svalue,Sname) updates the values of existing scenario
%  parameters in the dataset for all QOIs. The values and names of the scenario parameters 
%  are specified by the input Svalue and Sname, respectively.
%
%  changeScenarioParameter(obj,Svalue,Sname，QOIindex) updates the values of existing scenario
%  parameters in the dataset for QOIs specified by QOIindex. The values and names of
%  the scenario parameters are specified by the input Svalue and Sname, respectively.

%  obj - A B2BDC.B2Bdataset.Dataset object
%  Svalue - A numerical array, each column correponds to one scenario parameter
%  Sname - A cell array, specifies which scenario parameters to be changed

if nargin < 4
   QOIindex = 1:obj.Length;
elseif max(QOIindex > obj.Length)
   error('Wrong QOI index');
end
if length(QOIindex) ~= size(Svalue,1)
   error('Wrong dimension of scenario parameter value input');
end
[~,oldName] = obj.getScenarioParameter;
[~,id] = intersect(Sname,oldName,'stable');
if length(id) ~= length(Sname)
   error('The input scenario parameter is not found')
end
dsUnits = obj.DatasetUnits.Values;
for i = 1:length(QOIindex)
   tmpsv = dsUnits(QOIindex(i)).ScenarioParameter;
   tmpsv.Value(id) = Svalue(i,:);
   oldUnit = dsUnits(QOIindex(i));
   dsName = oldUnit.Name;
   newUnit = B2BDC.B2Bdataset.DatasetUnit(oldUnit.Name,oldUnit.SurrogateModel,oldUnit.LowerBound,...
      oldUnit.UpperBound,oldUnit.ObservedValue,tmpsv);
   obj.DatasetUnits = obj.DatasetUnits.replace('Name',dsName,newUnit);
%    dsUnits(QOIindex(i)).ScenarioParameter.Value(id) = Svalue(i,:);
end
% ds = generateDataset(obj.Name);
% for i = 1:numel(dsUnits)
%    ds.addDSunit(dsUnits(i));
% end