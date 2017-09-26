function [RMSE]= loocv(y,X)
%--------------------------------------------------------------------------------
% PURPOSE: Performs Leave One Out Cross-Validation on a linear model
%--------------------------------------------------------------------------------
% INPUTS: y = N-by-1 dependent variable
%         X = N-by-k independent variables
%--------------------------------------------------------------------------------
% OUTPUT: RMSE = Root Mean Squared Error (N-by-1)
  
diff = NaN(length(y),1);
for i = 1:length(y)
  Xlo      = X;
  Xlo(i,:) = [];
  ylo      = y;
  ylo(i)   = [];
  
  betas = (Xlo'*Xlo)\(Xlo'*ylo);
  yhat  = X*betas;
  diff(i) = yhat(i)-y(i);
end

RMSE = sqrt(sum(diff.^2)/length(y));
return
