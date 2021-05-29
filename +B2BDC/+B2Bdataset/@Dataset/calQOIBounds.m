function qoiPos = calQOIBounds(obj, index, opt)
%  qoiPos = calQOIBounds(obj) calculates the posterior uncertainty bound of all QOIs
%  in the dataset. The option variable is set to be the default B2BDC option
%
%  qoiPos = calQOIBounds(obj, index) calculates the posterior uncertainty bound of the
%  QOIs specified by index. The option variable is set to be the default B2BDC option.
%   
%  qoiPos = calQOIBounds(obj, index, opt) calculates the posterior uncertainty bound of
%  the QOIs specified by index. The option variable is specified by opt.

%  Created: Nov 30, 2015     Wenyu Li
%  Modified: Nov 16, 2016     Wenyu Li

if nargin < 2
   index = [];
end
if nargin < 3
   opt = generateOpt;
   opt.Display = true;
   opt.ExtraLinFraction = -1;
end
pflag = opt.Prediction;
nUnit = obj.Length;
if index > nUnit
    error('Index cannot exceed the number of QOIs in the dataset.');
end
allName = {obj.DatasetUnits.Values.Name};
if isempty(index)
   id = 1:nUnit;
elseif iscell(index)
   [~,~,id] = intersect(index, allName, 'stable');
elseif isnumeric(index)
   id = index;
else
   error('Wrong input index format')
end
selectName = allName(id);
qbd1 = zeros(length(id),2);
qbd2 = zeros(length(id),2);
% if opt.Display
   h = waitbar(0,'Calculating the Posterior QOI Bounds');
% end
for i = 1:length(id)
   tep_QOI = obj.DatasetUnits.Values(id(i)).SurrogateModel;
   QOIrange = obj.predictQOI(tep_QOI, opt);
   if strcmp(pflag,'both')
      qbd1(i,1) = QOIrange.min(1);
      qbd1(i,2) = QOIrange.max(2);
      qbd2(i,1) = QOIrange.min(2);
      qbd2(i,2) = QOIrange.max(1);
   else
      qbd1(i,1) = QOIrange.min;
      qbd1(i,2) = QOIrange.max;
   end
%    if opt.Display
      waitbar(i/length(id),h);
%    end
end
% if opt.Display
   close(h);
% end
qoiPos.varName = selectName;
switch pflag
   case 'both'
      qoiPos.OuterBound = qbd1;
      qoiPos.InnerBound = qbd2;
   case 'inner'
      qoiPos.InnerBound = qbd1;
   case 'outer'
      qoiPos.OuterBound = qbd1;
end

if opt.Display
   figure('NumberTitle','off',...
      'Units','normalized','Position',[0.1,0.1,0.8,0.8],...
      'Name','Posterior QoI bounds');
   bds = obj.calBound;
   obs = obj.calObserve;
   dub = bds(:,2)-obs;
   dlb = obs-bds(:,1);
   errorbar(id,obs(id),dlb(id),dub(id),'ok','MarkerSize',1e-5,'LineWidth',2);
   hold on
   switch pflag
      case 'inner'
         mbd = mean(qbd1,2);
         dx = qbd1(:,2)-mbd;
         errorbar(id-0.2,mbd,dx,'ob','MarkerSize',1e-5,'LineWidth',2);
         legend('Original QoI bounds','B2B predicted inner bounds')
      case 'outer'
         mbd = mean(qbd1,2);
         dx = qbd1(:,2)-mbd;
         errorbar(id+0.2,mbd,dx,'or','MarkerSize',1e-5,'LineWidth',2);
         legend('Original QoI bounds','B2B predicted outer bounds')
      otherwise
         mbd1 = mean(qbd1,2);
         dx1 = qbd1(:,2)-mbd1;
         mbd2 = mean(qbd2,2);
         dx2 = qbd2(:,2)-mbd2;
         errorbar(id-0.2,mbd2,dx2,'ob','MarkerSize',1e-5,'LineWidth',2);
         errorbar(id+0.2,mbd1,dx1,'or','MarkerSize',1e-5,'LineWidth',2);
         legend('Original QoI bounds','B2B predicted inner bounds','B2B predicted outer bounds')
   end
   hold off
   xlabel('QoI index')
   ylabel('Uncertainty bounds')
   set(gca,'XLim',[0,max(id)+1],'FontSize',17,'FontWeight','Bold','Position',[0.08,0.15,0.85,0.75]);
end