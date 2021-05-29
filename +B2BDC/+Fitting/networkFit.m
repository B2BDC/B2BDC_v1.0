function y = networkFit(X,Y,vars,HL,nCV)
% Fit a network model of the given data X(nSample-by-nVariable) and
% Y(nSample-by-1) with variable information in the
% B2BDC.B2Bvariables.VariableList object vars. 

% output:
% A B2BDC.B2Bmodels.NNModel object

if isempty(X)
   net = feedforwardnet(HL);
   net.trainFcn = 'trainbr';
   nVar = vars.Length;
   net = configure(net,randn(nVar,2),randn(1,2));
   while all(size(net.IW{1}) ~= [HL(1) nVar])
      net = configure(net,randn(nVar,2),randn(1,2));
   end
   y = B2BDC.B2Bmodels.NNModel(net,vars);
   return;
end
netEpo = 2e3;
netGrad = 1e-6;
netReg = 0;
nSample = size(X,1);
if nCV == 1
   net = feedforwardnet(HL);
   net.trainFcn = 'trainbr';
   net.performParam.regularization = netReg;
   net.trainParam.epochs = netEpo;
   net.trainParam.min_grad = netGrad;
%    net = configure(net,randn(nVar,2),randn(1,2));
   net.divideParam.trainRatio = 70/100;
   net.divideParam.valRatio = 15/100;
   net.divideParam.testRatio = 15/100;
   net.trainParam.showWindow = false;
   [net,tr] = train(net,X',Y');
   dy = abs(net(X')-Y');
   err.absMax = max(dy);
   err.absAvg = mean(dy);
   dy = dy./abs(Y');
   err.relMax = max(dy);
   err.relAvg = mean(dy);
   y = B2BDC.B2Bmodels.NNModel(net,vars,err);
else
   nTest = floor(nSample/nCV);
   err.absMax = 0;
   err.absAvg = 0;
   err.relMax = 0;
   err.relAvg = 0;
   for i = 1:nCV
      if i~=nCV
         idTest = (1:nTest)+(i-1)*nTest;
      else
         idTest = (nCV-1)*nTest+1:nSample;
      end
      idTrain = setdiff(1:nSample,idTest);
      xTrain = X(idTrain,:);
      yTrain = Y(idTrain);
      xTest = X(idTest,:);
      yTest = Y(idTest);
      net = feedforwardnet(HL);
      net.trainFcn = 'trainbr';
      net.performParam.regularization = netReg;
      net.trainParam.epochs = netEpo;
      net.trainParam.min_grad = netGrad;
      net.divideParam.trainRatio = 70/100;
      net.divideParam.valRatio = 15/100;
      net.divideParam.testRatio = 15/100;
      net.trainParam.showWindow = false;
      [net,tr] = train(net,xTrain',yTrain');
      dy = abs(net(xTest')-yTest');
      err.absMax = max(err.absMax,max(dy));
      err.absAvg = err.absAvg+mean(dy);
      dy = dy./abs(yTest');
      err.relMax = max(err.relMax,max(dy));
      err.relAvg = err.relAvg+mean(dy);
   end
   net = feedforwardnet(HL);
   net.trainFcn = 'trainbr';
   net.performParam.regularization = netReg;
   net.trainParam.epochs = netEpo;
   net.trainParam.min_grad = netGrad;
   net.divideParam.trainRatio = 70/100;
   net.divideParam.valRatio = 15/100;
   net.divideParam.testRatio = 15/100;
   net.trainParam.showWindow = false;
   [net,tr] = train(net,X',Y');
   y = B2BDC.B2Bmodels.NNModel(net,vars,err);
end
