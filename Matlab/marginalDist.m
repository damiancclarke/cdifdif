function [Xdist,ndist,maxdist] = marginalDist(y,X,dist,step,tlimit)
%--------------------------------------------------------------------------------
% PURPOSE: Tests up until marginal spillover is insignificant
%--------------------------------------------------------------------------------
% INPUTS:  y    = N-by-1 dependent variable
%          X    = N-by-k baseline independent variables
%          dist = N-by-1 distance to treatment
%          step = scalar distance step size
%--------------------------------------------------------------------------------
% OUTPUTS: Xdist   = N-by-(k+ndist) independent variables with new distance dummies
%          ndist   = Number of distance variables
%          maxdist = Maximum spillover distance
%--------------------------------------------------------------------------------

  tcrit   = 100;
  dl      = 0;
  [N,var] = size(X);
  Xdist   = X;
  margin  = var+1;
  ndist   = 0;
  
  while tcrit > tlimit
    dln   = dl+step;
    dnew  = dist>dl&dist<=dln;
    Xdist = [Xdist dnew];

    betas = (Xdist'*Xdist)\(Xdist'*y);
    uhat  = Xdist*betas-y;
    ses   = diag(sqrt((uhat'*uhat)/(N-length(betas))*inv(Xdist'*Xdist)));
    ts    = betas./ses;
    
    dl = dl+step;
    tcrit   = abs(ts(margin,1));
    margin  = margin+1;
    ndist   = ndist+1;
  end

  [n1,n2]=size(Xdist);
  Xdist=Xdist(:,1:n2-1);
  ndist=ndist-1;
  maxdist=dl-step;
