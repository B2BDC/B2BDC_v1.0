function err = showFittingError(obj,errorType)

% show distribution of fitting errors

%  Created: June 16, 2019     Wenyu Li

% errorType (not case sensitive): abs/absolute - absolute error
%                                 rel/relative - relative error


err = zeros(obj.Length,2);
if nargin < 2
   errorType = 'abs';
end
for i = 1:obj.Length
   errStat = obj.DatasetUnits.Values(i).SurrogateModel.ErrorStats;
   if strcmpi(errorType,'abs') || strcmpi(errorType,'absolute')
      err(i,:) = [errStat.absMax errStat.absAvg];
      titles = {'Maximum infinity-norm error','Average infinity-norm error'};
   elseif strcmpi(errorType,'rel') || strcmpi(errorType,'relative')
      err(i,:) = [errStat.relMax errStat.relAvg];
      titles = {'Maximum 2-norm error','Average 2-norm error'};
   else
      error('Wrong input error type');
   end
end
f = figure('Name','Error distribution');
f.Units = 'normalized';
f.Position = [0.13 0.09 0.52 0.77];
ww = 0.86;
hh = 0.38;
bp = [0.94-hh 0.06];
lp = 0.07;
for i = 1:2
   subplot(2,1,i);
   histogram(err(:,i));
   xx = get(gca,'XLim');
   xx(1) = 0;
   set(gca,'Position',[lp bp(i) ww hh],'LineWidth',1.5,'FontSize',15,'YTick',[],'XLim',xx);
   hold on
   plot(mean(err(:,i))*[1 1],get(gca,'YLim'),'k--','LineWidth',2);
   hold off
   text(0.3,1.08,titles{i},'Units','normalized','FontSize',20)
   tt = ["max: "+num2str(max(err(:,i)),'%.4f'); "min: "+num2str(min(err(:,i)),'%.4f'); "avg: "+num2str(mean(err(:,i)),'%.4f')];
   text(0.86,0.85,tt,'Units','normalized','FontSize',15);
end
