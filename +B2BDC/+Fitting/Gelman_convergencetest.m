function [R,Rx,I,det_V,det_W] = Gelman_convergencetest(x,showPlot)

% Calculate Gelman inference from multiple samples chains of a vector x.

% Modified: June 28, 2019     Wenyu Li

if nargin < 2
   showPlot = false;
end
frac = 0.1:0.1:1;
nPlot = length(frac);
det_V = zeros(nPlot,1);
det_W = zeros(nPlot,1);
I = zeros(nPlot,1);
nTest = round(size(x,1)*frac);
for i = 1:nPlot
   xx = x(1:nTest(i),:,:);
   [n,nv,m] = size(xx);
   W = zeros(nv);
   for j = 1:m
      W = W+cov(xx(:,:,j));
   end
   W = W/m;
   mx = mean(xx);
   mx = permute(mx,[3 2 1]);
   B = cov(mx);
   V = (n-1)/n*W + (m+1)/m*B;
   det_V(i) = det(V);
   det_W(i) = det(W);
   lambda = max(eig(W\B));
   I(i) = (n-1)/n+(m+1)/m*lambda;
end
R = sqrt(I(end));
if showPlot
   f = gcf;
   if strcmp(f.Units,'pixels')
      set(f,'defaultAxesColorOrder',[0 0 0; 1 0 0]);
      clf;
      f.Position = [270 320 1250 580];
   end
   xx = m*nTest;
   yyaxis left
   hold on
   plot(xx,sqrt(I),'-k','LineWidth',2.5);
   plot([0 m*1.1*nTest(end)],[1 1],'--k','LineWidth',2.5);
   hold off
   limsy = get(gca,'YLim');
   set(gca,'YLim',[0 limsy(2)]);
   yyaxis right
   hold on
   hlg(1) = plot(xx,log10(det_V),'b-','LineWidth',2.5);
%    plot(xx,sqrt(est_sigma),'r-','LineWidth',2.5);
   hlg(2) = plot(xx,log10(det_W),'g-','LineWidth',2.5);
   hold off
%    lg = legend('Potential reduction','V hat','sigma','W');
%    lg.Position = [0.8 0.2 0.1368 0.1405];
   set(gca,'Position',[0.08 0.1 0.8,0.85],'FontSize',12,'LineWidth',1.5,...
      'XLim',[0,m*1.1*nTest(end)],'XTick',m*nTest,'XScale','linear');
   box on;
   xlabel('Number of sample');
   lg = legend(hlg,'log |V|','log |W|','Location','best');
   lg.FontSize = 18;
   lg.Box = 'off';
   text(-0.062,0.4351,'MPSRF','Units','normalized','FontSize',18,'Rotation',90);
   text(1.084,0.2363,'Logarithmic determinant','Units','normalized','FontSize',18,'Rotation',90);
end


[n,nv,m] = size(x);
Rx = zeros(nv,1);
for i = 1:nv
   xx = permute(x(:,i,:),[1 3 2]);
   mx = mean(xx);
   W = var(xx);
   B = var(mx);
   sig = (n-1)*n*W+B;
   Rx(i) = (m+1)/m*sig/W-(n-1)/m/n;
end