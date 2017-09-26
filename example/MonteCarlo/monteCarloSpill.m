% monteCarloSpill.m             damiancclarke              yyyy-mm-dd:2017-01-01
%---|----1----|----2----|----3----|----4----|----5----|----6----|----7----|----8
%--
%--  Simulate estimation of spillovers where spillovers are captured by the algo
%-- rithm described in "Estimating Difference-in-Differences in the Presence of
%-- Spillovers".  This tests three models using a constant bandwidth and Leave-
%-- One-Out Cross-Validation.  It also checks models with stratified k-fold
%-- Cross-Validation.

diary on
clear
rng(272)

sims    = 2500;
N       = 1000;
effects = [3 10 5 4 3 2];
  
treat = zeros(N,1);
tstart=0.8*N+1
treat(tstart:N)=1;

output   = NaN(12,12);
outputkf = NaN(12,12);

dd=0;
dfactors = [2.5*N/1000 5*N/1000 12.5*N/1000];

for dsim = dfactors
  dd = dd + 1;
  dist = transpose(([1:N])/dsim);
  dist = dist.*(1-treat);
  
  d20  = dist>0&dist<=20;
  d10  = dist>0&dist<=10;

  tvalsB   = NaN(sims,9);
  bEst     = NaN(sims,9);
  hvals    = NaN(sims,9);
  tvalsBkf = NaN(sims,9);
  bEstkf   = NaN(sims,9);
  hvalskf  = NaN(sims,9);

  fprintf('Simulation distance %d...\n',dsim)
  %-------------------------------------------------------------------------------
  %--- (1) Spillovers
  %-------------------------------------------------------------------------------
  for round = 1:9
    error = mod(round,3);
    if error==0
      error=5;
    end

    RMSEave        = zeros(25,2);
    RMSEave(:,1)   = 1:25;
    RMSEavekf      = zeros(25,2);
    RMSEavekf(:,1) = 1:25;
    hold on
    
    for simul = 1:sims
      if round>=1&round<=3
        %MODEL 1
        dmodel  = [0 5 10 15 20];
	loc = round+1;
      elseif round>=4&round<=6
        %MODEL 2
        dmodel  = [0 2 9 17 20];
	loc = round+2;
      else
        %MODEL 3
        y = effects(1) + effects(2)*treat + 10*exp(-dist/5).*d10 + error*randn(N,1);
        ceff = sum(10*exp(-dist(1:50)/5))/length(dist(1:50));
	loc = round+3;
      end
  
      if round<7
        effects = [3 10 5 4 3 2];
        d1    = dist>dmodel(1)&dist<=dmodel(2);
        d2    = dist>dmodel(2)&dist<=dmodel(3);
        d3    = dist>dmodel(3)&dist<=dmodel(4);
        d4    = dist>dmodel(4)&dist<=dmodel(5);
        y = effects(1) + effects(2)*treat + effects(3)*d1 + ...
  	  effects(4)*d2 + effects(5)*d3 + effects(6)*d4 + error*randn(N,1);
        ceff  = ((dmodel(2)-dmodel(1))/dmodel(5))*effects(3)+...
    	      ((dmodel(3)-dmodel(2))/dmodel(5))*effects(4)+...
    	      ((dmodel(4)-dmodel(3))/dmodel(5))*effects(5)+...
    	      ((dmodel(5)-dmodel(4))/dmodel(5))*effects(6);
      end
      
      Xprop = [ones(N,1) treat];
      
      j          = 1;
      RMSEmin    = Inf;
      RMSEminkf  = Inf;
      optdist    = 0;
      optdistkf  = 0;
      Xest       = Xprop;
      Xestkf     = Xprop;
  
      RMSEres = NaN(25,2);
      for i = linspace(1,25,25)
        [Xspill,nspill,mdist] = marginalDist(y,Xprop,dist,i);
        RMSE                = loocv(y,Xspill);
	nvars = ndims(Xspill);
	if nvars>2
	  strat=Xspill(:,3);
	else
	  strat=Xspill(:,2);
        end
	regf=@(XTRAIN,ytrain,XTEST)(XTEST*regress(ytrain,XTRAIN));
        RMSEkfold = crossval('mse',Xspill,y,'predfun',regf,'stratify',strat);

        if RMSE<RMSEmin
  	  RMSEmin = RMSE;
  	  optdist = i;
  	  Xest    = Xspill;
  	  maxdist = mdist;
        end
        if RMSEkfold<RMSEminkf
  	  RMSEminkf = RMSEkfold;
  	  optdistkf = i;
  	  Xestkf    = Xspill;
  	  maxdistkf = mdist;
        end
        RMSEres(j,1)  = i;
        RMSEres(j,2)  = RMSE;
        RMSEave(j,2)  = RMSEave(j,2)   + RMSE/sims;
        RMSEavekf(j,2)= RMSEavekf(j,2) + RMSEkfold/sims;
        j = j + 1;
      end

      h1=plot(RMSEres(:,1),RMSEres(:,2),'LineWidth',1,'color',[0.8 0.8 0.8]);

      %%LEAVE-ONE-OUT CV VERSION
      betasFinal = (Xest'*Xest)\Xest'*y;
      uhat       = Xest*betasFinal-y;
      sesFinal   = diag(sqrt((uhat'*uhat)/(N-length(betasFinal))*inv(Xest'*Xest)));  
      sdist      = dist>0&dist<=maxdist;
      Xsimp      = [ones(N,1) treat sdist];
      betasSimp  = (Xsimp'*Xsimp)\Xsimp'*y;
      uhat       = (Xsimp*betasSimp)-y;
      sesSimp    = diag(sqrt((uhat'*uhat)/(N-length(betasSimp))*inv(Xsimp'*Xsimp)));
      
      tB = (betasFinal(2)-effects(2))/sesFinal(2);
      tC = (betasSimp(3)-ceff)/sesSimp(3);
      
      tvalsB(simul,round)=tB;
      bEst(simul,round)   =betasFinal(2);
      hvals(simul,round)  =optdist;

      %%K-FOLD CV VERSION
      betasFinalkf = (Xestkf'*Xestkf)\Xestkf'*y;
      uhat         = Xestkf*betasFinalkf-y;
      sesFinalkf   = diag(sqrt((uhat'*uhat)/(N-length(betasFinalkf))*inv(Xestkf'*Xestkf)));  
      sdistkf      = dist>0&dist<=maxdistkf;
      Xsimpkf      = [ones(N,1) treat sdistkf];
      betasSimpkf  = (Xsimpkf'*Xsimpkf)\Xsimpkf'*y;
      uhat         = (Xsimpkf*betasSimpkf)-y;
      sesSimpkf    = diag(sqrt((uhat'*uhat)/(N-length(betasSimpkf))*inv(Xsimpkf'*Xsimpkf)));
      
      tBkf = (betasFinalkf(2)-effects(2))/sesFinalkf(2);
      tCkf = (betasSimpkf(3)-ceff)/sesSimpkf(3);
      
      tvalsBkf(simul,round)=tBkf;
      bEstkf(simul,round)   =betasFinalkf(2);
      hvalskf(simul,round)  =optdistkf;

		    
      count = (round-1)*sims+simul;
      fprintf('%d...', count)
      fprintf('Simulation %d, round %d, Spillover %5.3f \n',simul,round,maxdist)
    end
    output(loc,dd*4-3) = mean(hvals(:,round));
    output(loc,dd*4-2) = mean(bEst(:,round));
    output(loc,dd*4-1) = std(bEst(:,round));
    output(loc,dd*4+0) = mean(abs(tvalsB(:,round))>1.96);

    outputkf(loc,dd*4-3) = mean(hvalskf(:,round));
    outputkf(loc,dd*4-2) = mean(bEstkf(:,round));
    outputkf(loc,dd*4-1) = std(bEstkf(:,round));
    outputkf(loc,dd*4+0) = mean(abs(tvalsBkf(:,round))>1.96);
		    
    h2 = plot(RMSEave(:,1),RMSEave(:,2),'-o','LineWidth',3,'MarkerSize',8,'Color',[0.2 0.6 0.6]);
	 hXLabel = xlabel('Bandwidth ($h$)','Interpreter','latex','FontSize',12);
         hYLabel = ylabel('RMSE $CV(h)$'   ,'Interpreter','latex','FontSize',12);
    set(gca, ...
	'Box'         , 'off'     , ...
	'TickDir'     , 'out'     , ...
	'TickLength'  , [.02 .02] , ...
	'XMinorTick'  , 'on'      , ...
	'YMinorTick'  , 'on'      , ...
	'XColor'      , [.3 .3 .3], ...
	'YColor'      , [.3 .3 .3], ...
	'XTick'       , 1:2:30    , ...
	'LineWidth'   , 1         );
    hleg = legend([h1 h2],{'Single Simulation','Average'},'Location','NorthEast');
    [minRMSE,argminRMSE] = min(RMSEave(:,2));
    str1 = sprintf('$$\\min_{ h}$$ RMSE $CV(h) =$ %3.5g',minRMSE);
    str2 = sprintf('arg$$\\min_{ h}$$ RMSE $CV(h) =$ %3.2g',argminRMSE);
    str  =  {str1,str2}
    DataX = interp1( [0 1], xlim(), 0.1);
    DataY = interp1( [0 1], ylim(), 0.9);
    text(DataX,DataY,str,'Interpreter','latex','FontSize',12);
    fname = strcat('results/RMSEplot_r',num2str(round),'_sim',num2str(dd))
    print(fname,'-depsc')
    hold off
    clf
		      
    plot(RMSEave(:,1),RMSEave(:,2),'-o','LineWidth',2,'Color',[0.2 0.6 0.6]);
	 hXLabel = xlabel('Bandwidth ($h$)','Interpreter','latex','FontSize',12);
         hYLabel = ylabel('RMSE $CV(h)$'   ,'Interpreter','latex','FontSize',12);
    hold on
    plot(RMSEavekf(:,1),RMSEavekf(:,2),'-s','LineWidth',2,'Color',[1 0.2 0.2]);
    legend('Leave-One-Out Cross Validation','k-fold Cross-Validation',...
	   'Location','NorthEast');
    set(gca, ...
	'Box'         , 'off'     , ...
	'TickDir'     , 'out'     , ...
	'TickLength'  , [.02 .02] , ...
	'XMinorTick'  , 'on'      , ...
	'YMinorTick'  , 'on'      , ...
	'XColor'      , [.3 .3 .3], ...
	'YColor'      , [.3 .3 .3], ...
	'XTick'       , 1:2:30    , ...
	'LineWidth'   , 1         );
    [minRMSE,argminRMSE] = min(RMSEave(:,2));
    [minRMSEkf,argminRMSEkf] = min(RMSEavekf(:,2));
    str1 = sprintf('arg$$\\min_{ h}$$ RMSE Leave-One-Out $CV(h) =$ %3.2g',argminRMSE);
    str2 = sprintf('arg$$\\min_{ h}$$ RMSE k-fold $CV(h) =$ %3.2g',argminRMSEkf);
    str  =  {str1,str2}
    DataX = interp1( [0 1], xlim(), 0.07);
    DataY = interp1( [0 1], ylim(), 0.77);
    text(DataX,DataY,str,'Interpreter','latex','FontSize',12);
    fname = strcat('results/RMSEmean_r',num2str(round),'_sim',num2str(dd))
    print(fname,'-depsc')
    hold off
    clf
  end

  %-------------------------------------------------------------------------------
  %--- (2) Naive models
  %-------------------------------------------------------------------------------
  betasnaive=NaN(sims,3);
  tvalsnaive=NaN(sims,3);
  error = 3;
  
  for round = 1:3
    if round==1
      %MODEL 1
      dmodel  = [0 5 10 15 20];
    elseif round==2
      %MODEL 2
      dmodel  = [0 2 9 17 20];
    else
      %MODEL 3
      y = repmat(effects(1) + effects(2)*treat + 10*exp(-dist/5).*d10,1,sims) + ...
  	       error*randn(N,sims);
    end
  
    if round<3
      effects = [3 10 5 4 3 2];
      d1    = dist>dmodel(1)&dist<=dmodel(2);
      d2    = dist>dmodel(2)&dist<=dmodel(3);
      d3    = dist>dmodel(3)&dist<=dmodel(4);
      d4    = dist>dmodel(4)&dist<=dmodel(5);
      y = repmat(effects(1) + effects(2)*treat + effects(3)*d1 + ...
  	       effects(4)*d2 + effects(5)*d3 + effects(6)*d4,1,sims) + ...
  	       error*randn(N,sims);
      end
  
  
    Xprop = [ones(N,1) treat];
    betas = (Xprop'*Xprop)\Xprop'*y;
    uhat  = Xprop*betas-y;
    ses   = NaN(2,sims);
    for i=1:sims
      ses(:,i) = diag(sqrt((uhat(:,i)'*uhat(:,i))/(N-2)*inv(Xprop'*Xprop)));
    end
    betasnaive(:,round) = betas(2,:)';
    tvalsnaive(:,round) = abs((betas(2,:)-10)./ses(2,:))';
    loc = round*3-2+(round-1);
    output(loc,dd*4-2) = mean(betasnaive(:,round));
    output(loc,dd*4-1) = std(betasnaive(:,round));
    output(loc,dd*4+0) = mean(tvalsnaive(:,round)>1.96);
    outputkf(loc,dd*4-2) = mean(betasnaive(:,round));
    outputkf(loc,dd*4-1) = std(betasnaive(:,round));
    outputkf(loc,dd*4+0) = mean(tvalsnaive(:,round)>1.96);
  end
end

csep = {'&';'&';'&';'&';'&';'&';'&';'&';'&';'&';'&';'&'};
lsep = {'\\';'\\';'\\';'\\';'\\';'\\';'\\';'\\';'\\';'\\';'\\';'\\'};
nams = {'\textbf{Model 1}&&&&&&&&&&&\\Naive&';...
	'$\varepsilon=1$&';'$\varepsilon=2$&';'$\varepsilon=5$&';...
	'\textbf{Model 2}&&&&&&&&&&&\\Naive&';...
	'$\varepsilon=1$&';'$\varepsilon=2$&';'$\varepsilon=5$&';...
	'\textbf{Model 3}&&&&&&&&&&&\\Naive&';...
	'$\varepsilon=1$&';'$\varepsilon=2$&';'$\varepsilon=5$&'};

tout = dataset(nams,num2str(output(:,1),'%5.3f'),csep,...
	       num2str(output(:,2),'%5.3f'),csep,...
	       num2str(output(:,3),'%5.3f'),csep,...
	       num2str(output(:,4),'%5.3f'),csep,...
	       num2str(output(:,5),'%5.3f'),csep,...
	       num2str(output(:,6),'%5.3f'),csep,...
	       num2str(output(:,7),'%5.3f'),csep,...
	       num2str(output(:,8),'%5.3f'),csep,...
	       num2str(output(:,9),'%5.3f'),csep,...
	       num2str(output(:,10),'%5.3f'),csep,...
	       num2str(output(:,11),'%5.3f'),csep,...
	       num2str(output(:,12),'%5.3f'),lsep)
export(tout,'file','results/MCbetaVals.tex','Delimiter',' ',...
       'WriteVarNames',false)
	     
tout = dataset(nams,num2str(outputkf(:,1),'%5.3f'),csep,...
	       num2str(outputkf(:,2),'%5.3f'),csep,...
	       num2str(outputkf(:,3),'%5.3f'),csep,...
	       num2str(outputkf(:,4),'%5.3f'),csep,...
	       num2str(outputkf(:,5),'%5.3f'),csep,...
	       num2str(outputkf(:,6),'%5.3f'),csep,...
	       num2str(outputkf(:,7),'%5.3f'),csep,...
	       num2str(outputkf(:,8),'%5.3f'),csep,...
	       num2str(outputkf(:,9),'%5.3f'),csep,...
	       num2str(outputkf(:,10),'%5.3f'),csep,...
	       num2str(outputkf(:,11),'%5.3f'),csep,...
	       num2str(outputkf(:,12),'%5.3f'),lsep)
export(tout,'file','results/MCbeta-kfold.tex','Delimiter',' ',...
       'WriteVarNames',false)

diary off
