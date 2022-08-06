function [ypred,FVU,mdl] = fitlm_cv_saccadeResponse2(T,spec)

% This is for predicting the change in firing rate, and takes care of the
% nonlinearity

numK = 5;
cp = cvpartition(size(T,1),'KFold',numK);
e = 1;
ypred = zeros(size(T,1),1);
for i=1:numK
    exc = e:e+cp.TestSize(i)-1;
    e = e+cp.TestSize(i);
    mdl = fitlm(T,spec,'Intercept',false,'Exclude',exc);
    y = predict(mdl,T(exc,:));
    ypred(exc) = y;
end

% Based on the baseline FR, figure out the maximum allowed decrease in
% firing rate
maxDecrease = -T.baselineGrating;
% if ypred predicts decrease greater than maxDecrease, set to maxDecrease
ypred(ypred<maxDecrease)=maxDecrease(ypred<maxDecrease);

SSerr = sum((T.gratingChange-ypred).^2);
y_mean = mean(T.gratingChange);
SStot = sum((T.gratingChange-y_mean).^2);
FVU = SSerr/SStot;
mdl = fitlm(T,spec,'Intercept',false);