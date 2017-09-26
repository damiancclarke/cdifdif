function [RMSE]= kfoldcv(y,X,k)
%--------------------------------------------------------------------------------
% PURPOSE: Performs k-fold Cross-Validation on a linear model
%--------------------------------------------------------------------------------
% INPUTS: y = N-by-1 dependent variable
%         X = N-by-k independent variables
%         k = Number of folds
%--------------------------------------------------------------------------------
% OUTPUT: RMSE = Root Mean Squared Error (N-by-1)
  
diff = NaN(length(y),1);
rnum = ceil(rand(length(y),1)*k);

for i = 1:k
  Xlo            = X;
  Xlo(rnum~=i,:) = [];
  ylo            = y;
  ylo(rnum~=i)   = [];
  
  betas = (Xlo'*Xlo)\(Xlo'*ylo);
  yhat  = X*betas;
  diff(rnum==i) = yhat(rnum==i)-y(rnum==i);
end

RMSE = sqrt(sum(diff.^2));
return
