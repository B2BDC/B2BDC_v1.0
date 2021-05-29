function  y = eval(obj,x, varObj)
   % y = eval(obj, x) evaluates the model output at input design points specified by X.
   % Each row of X corresponds to a design point.

   % y = eval(obj, x, varObj) evaluates the model output at input design points specified
   % by X and varObj. Each row of X corresponds to a design point. The column of X is
   % specified by the given input varObj

% Created: July 21, 2015    Wenyu Li
unitModel = obj.SurrogateModel;
if nargin > 2
   y = unitModel.eval(x, varObj);
else
   y = unitModel.eval(x);
end
