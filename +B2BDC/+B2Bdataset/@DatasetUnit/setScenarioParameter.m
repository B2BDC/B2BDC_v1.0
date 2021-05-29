function newOBJ = setScenarioParameter(oldOBJ,sv,sname)
   % newOBJ = setScenarioParameter(oldOBJ, sv, sname) sets the scenario parameter of
   % the dataset unit. The values and names of the scenario parameters are specified by sv
   % and sname, respectively.

% Created: Sep 9, 2018     Wenyu Li
if ~iscell(sname)
   sname = {sname};
end
if length(sv) ~= length(sname)
   error('Wrong input dimension');
else
   newOBJ = oldOBJ;
   newOBJ.ScenarioParameter.Value = sv;
   newOBJ.ScenarioParameter.Name = sname; 
end