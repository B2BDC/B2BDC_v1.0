function clearScenarioParameter(obj)
%  clearScenarioParameter(obj) clears all existing scenario parameters in the dataset.

%  obj - A B2BDC.B2Bdataset.Dataset object
%  Svalue - A numerical array, each column correponds to one scenario parameter
%  Sname - A cell array, specifies which scenario parameters to be changed

dsUnits = obj.DatasetUnits.Values;
for i = 1:obj.Length
   tmpsv.Value = [];
   tmpsv.Name = [];
   oldUnit = dsUnits(i);
   dsName = dsUnits(i).Name;
   newUnit = B2BDC.B2Bdataset.DatasetUnit(oldUnit.Name,oldUnit.SurrogateModel,oldUnit.LowerBound,...
      oldUnit.UpperBound,oldUnit.ObservedValue,tmpsv);
   obj.DatasetUnits = obj.DatasetUnits.replace('Name',dsName,newUnit);
end