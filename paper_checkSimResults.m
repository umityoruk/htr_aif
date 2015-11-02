load('simResults_s0p2.mat')
%load('simResults.mat')
load('truePSC.mat')

paper_loadParameters;
% GFR

GFRErrVS = simResults.GFRErrVS;
GFRErrHTR = simResults.GFRErrHTR;
GFRErrTRUE = simResults.GFRErrTRUE;

figure
imagesc(abs(GFRErrVS))
title('VS Error')

figure
imagesc(abs(simResults.GFRErrHTR))
title('HTR Error')


medErrVS = median(abs(GFRErrVS),1);
medErrHTR = median(abs(GFRErrHTR),1);
medErrTRUE = median(abs(GFRErrTRUE),1);

figure
plot(GFR(2:end),medErrVS(2:end),'b-')
hold on
plot(GFR(2:end),medErrHTR(2:end),'r-')
plot(GFR(2:end),medErrTRUE(2:end),'k--')
hold off
title('Median Error')
xlabel('GFR [ml/min]')
ylabel('% Error')
legend('VS','HTR','TRUE')


meanErrVS = mean(abs(GFRErrVS),1);
meanErrHTR = mean(abs(GFRErrHTR),1);
meanErrTRUE = mean(abs(GFRErrTRUE),1);

figure
plot(GFR(2:end),meanErrVS(2:end),'b-')
hold on
plot(GFR(2:end),meanErrHTR(2:end),'r-')
plot(GFR(2:end),meanErrTRUE(2:end),'k--')
hold off
title('Mean Error')
xlabel('GFR [ml/min]')
ylabel('% Error')
legend('VS','HTR','TRUE')


maxErrVS = max(abs(GFRErrVS),[],1);
maxErrHTR = max(abs(GFRErrHTR),[],1);
maxErrTRUE = max(abs(GFRErrTRUE),[],1);


figure
plot(GFR(2:end),maxErrVS(2:end),'b-')
hold on
plot(GFR(2:end),maxErrHTR(2:end),'r-')
plot(GFR(2:end),maxErrTRUE(2:end),'k--')
hold off
title('Max Error')
xlabel('GFR [ml/min]')
ylabel('% Error')
legend('VS','HTR','TRUE')



minErrVS = min(abs(GFRErrVS),[],1);
minErrHTR = min(abs(GFRErrHTR),[],1);
minErrTRUE = min(abs(GFRErrTRUE),[],1);

gStart = 2;
xvec = GFR(gStart:end);
xvec = [xvec,fliplr(xvec)];
yHTR = [minErrHTR(gStart:end),fliplr(maxErrHTR(gStart:end))];
yVS = [minErrVS(gStart:end),fliplr(maxErrVS(gStart:end))];
figure
h1 = fill(xvec,yHTR,'r');
hold on
h2 = fill(xvec,yVS,'b');
hold off
set(h1,'FaceAlpha',0.5);
set(h2,'FaceAlpha',0.5);



fitInfoMat = simResults.fitInfoMat;

trueAIF = truePSC.trueAIF;
trueCkidney = truePSC.Ckidney;
t_true = 0:length(trueAIF)-1;       


gfrInd = 8;
figure
rmse_htr_vec = zeros(size(fitInfoMat,1),length(GFR));
rmse_vs_vec = zeros(size(fitInfoMat,1),length(GFR));
for gg = 1:length(GFR)
for dd=1:size(fitInfoMat,1)
    fitInfo = fitInfoMat{dd,gg};
    t_htr = fitInfo.t_aortaHTR-1;
    t_vs = fitInfo.t_aortaVS-1;
    plot(0:length(trueAIF)-1,trueAIF,'k-','LineWidth',2)
    hold on
    
%     Cp_aortaHTR_smooth = smooth(fitInfo.Cp_aortaHTR,fitInfo.smoothBest,'loess');
%     plot(fitInfo.t_aortaHTR,Cp_aortaHTR_smooth,'r-','LineWidth',2)

    plot(t_htr,fitInfo.Cp_aortaHTR,'r-','LineWidth',2)
    plot(t_vs,fitInfo.Cp_aortaVS,'b-','LineWidth',2)
    rmse_htr_vec(dd,gg) = sqrt(sum((fitInfo.Cp_aortaHTR - interp1(t_true,trueAIF,t_htr)).^2)./length(t_htr));
    rmse_vs_vec(dd,gg) = sqrt(sum((fitInfo.Cp_aortaVS - interp1(t_true,trueAIF,t_vs)).^2)./length(t_vs));
%     hold off
%     pause
end
end
plot(t_true,trueAIF,'k-','LineWidth',2)
title('AIF estimation over 12 jitters')
xlabel('Time [s]')
ylabel('% Signal Change')
legend('TrueAIF','HTR','VS')

figure
plot(rmse_htr_vec(:,gfrInd),'r.')
hold on
plot(rmse_vs_vec(:,gfrInd),'b.')

fprintf(['Mean VS-AIF RMSE  = ' num2str(mean(mean(rmse_vs_vec,1),2)) '\n'])
fprintf(['Mean HTR-AIF RMSE = ' num2str(mean(mean(rmse_htr_vec,1),2)) '\n'])

gfrInd = 8;
figure
for dd=1:size(fitInfoMat,1)
    trueKidney = trueCkidney(:,gfrInd);
    fitInfo = fitInfoMat{dd,gfrInd};
    plot(0:length(trueKidney)-1,trueKidney,'k-','LineWidth',2)
    hold on
    %     xdata = [fitInfo.Cp_aortaHTR,fitInfo.t_aortaHTR];
    %     [CkidneyFit,Cart,Ctub] = ThreeCompartment(fitInfo.xFitHTR,xdata);
    %     plot(fitInfo.t_aortaHTR,CkidneyFit,'r-','LineWidth',2)
    
    
%     CkidneyHTR_smooth = smooth(fitInfo.CkidneyHTR,fitInfo.smoothBest,'loess');
%     plot(fitInfo.t_aortaHTR,CkidneyHTR_smooth,'r-','LineWidth',2)
    
    plot(fitInfo.t_aortaHTR-1,fitInfo.CkidneyHTR,'r-','LineWidth',2)
    plot(fitInfo.t_aortaVS-1,fitInfo.Ckidney,'b-','LineWidth',2)
%     hold off
%     pause
end
title('Ckidney estimation over 12 jitters')
xlabel('Time [s]')
ylabel('% Signal Change')
legend('TrueCkidney','HTR','VS')





gfrInd = 8;
figure
for dd=1:size(fitInfoMat,1)
    trueKidney = trueCkidney(:,gfrInd);
    fitInfo = fitInfoMat{dd,gfrInd};
    plot(0:length(trueKidney)-1,trueKidney,'k-','LineWidth',2)
    hold on
    xdata = [fitInfo.Cp_aortaHTR,fitInfo.t_aortaHTR];
    [CkidneyFit,Cart,Ctub] = ThreeCompartment(fitInfo.xFitHTR,xdata);
    plot(fitInfo.t_aortaHTR,CkidneyFit,'g.','LineWidth',2)
    xdata = [fitInfo.Cp_aortaVS,fitInfo.t_aortaVS];
    [CkidneyFit,Cart,Ctub] = ThreeCompartment(fitInfo.xFitVS,xdata);
    plot(fitInfo.t_aortaVS,CkidneyFit,'b.','LineWidth',2)
%     hold off
%     pause
end
title('Ckidney estimation over 12 jitters')
xlabel('Time [s]')
ylabel('% Signal Change')
legend('TrueCkidney','HTR','VS')

%%

resnormMatVS = simResults.resnormMatVS;
resnormMatHTR = simResults.resnormMatHTR;
resnormMatTRUE = simResults.resnormMatTRUE;

medResVS = median(abs(resnormMatVS),1);
medResHTR = median(abs(resnormMatHTR),1);
medResTRUE = median(abs(resnormMatTRUE),1);

figure
plot(GFR(2:end),medResVS(2:end),'b-')
hold on
plot(GFR(2:end),medResHTR(2:end),'r-')
plot(GFR(2:end),medResTRUE(2:end),'k--')
hold off
title('Median Resnorm')
xlabel('GFR [ml/min]')
ylabel('Resnorm')
legend('VS','HTR')


%% Scatter plot

GFRMatVS = simResults.GFRMatVS;
GFRMatHTR = simResults.GFRMatHTR;
GFRMatTRUE = simResults.GFRMatTRUE;

medGFRVS = median(abs(GFRMatVS),1);
medGFRHTR = median(abs(GFRMatHTR),1);
medGFRTRUE = median(abs(GFRMatTRUE),1);

% Least squares linear fit
A = [GFR(:) ones(length(GFR),1)];
B_VS = medGFRVS(:);
B_HTR = medGFRHTR(:);
B_TRUE = medGFRTRUE(:);
ab_VS = A\B_VS;
ab_HTR = A\B_HTR;
ab_TRUE = A\B_TRUE;
medLinearVS = A*ab_VS;
medLinearHTR = A*ab_HTR;
medLinearTRUE = A*ab_TRUE;

figure
plot(GFR,GFRMatVS(1,:),'bo','markerSize',10)
hold on
plot(GFR,GFRMatHTR(1,:),'rd','markerSize',10)

plot(GFR,GFR,'k--','LineWidth',2)
plot(GFR,GFRMatVS.','bo','markerSize',10)
plot(GFR,GFRMatHTR.','rd','markerSize',10)
plot(GFR,GFRMatTRUE.','ks','markerSize',10)
% plot(GFR,medGFRVS,'b-','LineWidth',2)
% plot(GFR,medGFRHTR,'r-','LineWidth',2)
plot(GFR,medLinearVS,'b-','LineWidth',2)
plot(GFR,medLinearHTR,'r-','LineWidth',2)
plot(GFR,medLinearTRUE,'k-','LineWidth',2)
hold off
legend('using VS-AIF','using HTR-AIF','Location','SouthEast')
xlabel('True GFR (ml/min)')
ylabel('Estimated GFR (ml/min)')
paper_setFigureProps;



%%
GFRStart = 2;
GFR_line = GFR(GFRStart):GFR(end);
medErrVS_line = interp1(GFR(GFRStart:end),medErrVS(GFRStart:end),GFR_line,'cubic');
medErrHTR_line = interp1(GFR(GFRStart:end),medErrHTR(GFRStart:end),GFR_line,'cubic');

figure
plot(GFR(GFRStart:end),medErrVS(GFRStart:end),'bo','markerSize',10,'LineWidth',2,'MarkerFaceColor','w')
hold on
plot(GFR(GFRStart:end),medErrHTR(GFRStart:end),'rd','markerSize',10,'LineWidth',2,'MarkerFaceColor','w')
plot(GFR(GFRStart:end),ones(size(GFR(GFRStart:end))).*5,'k--','LineWidth',2)
plot(GFR(GFRStart:end),ones(size(GFR(GFRStart:end))).*10,'k--','LineWidth',2)
plot(GFR(GFRStart:end),ones(size(GFR(GFRStart:end))).*15,'k--','LineWidth',2)
plot(GFR(GFRStart:end),ones(size(GFR(GFRStart:end))).*20,'k--','LineWidth',2)
plot(GFR(GFRStart:end),ones(size(GFR(GFRStart:end))).*25,'k--','LineWidth',2)
plot(GFR(GFRStart:end),ones(size(GFR(GFRStart:end))).*30,'k--','LineWidth',2)
plot(GFR(GFRStart:end),ones(size(GFR(GFRStart:end))).*35,'k--','LineWidth',2)
plot(GFR_line,medErrVS_line,'b-','LineWidth',2)
plot(GFR_line,medErrHTR_line,'r-','LineWidth',2)

plot(GFR(GFRStart:end),medErrVS(GFRStart:end),'bo','markerSize',10,'LineWidth',2,'MarkerFaceColor','w')
plot(GFR(GFRStart:end),medErrHTR(GFRStart:end),'rd','markerSize',10,'LineWidth',2,'MarkerFaceColor','w')

hold off
xlabel('GFR (ml/min)')
ylabel('% Error')
legend('using VS-AIF','using HTR-AIF','Location','NorthEast')
paper_setFigureProps;



%%
gfrInd = 8;
dd = 9;
figure
fitInfo = fitInfoMat{dd,gfrInd};
plot(0:length(trueAIF)-1,trueAIF*100,'k--','LineWidth',2)
hold on
plot(fitInfo.t_aortaVS-1,fitInfo.Cp_aortaVS*100,'bo-','LineWidth',2)
plot(fitInfo.t_aortaHTR-1,fitInfo.Cp_aortaHTR*100,'ro-','LineWidth',2)

% Cp_aortaHTR_smooth = smooth(fitInfo.Cp_aortaHTR,fitInfo.smoothBest,'loess');
% plot(fitInfo.t_aortaHTR=1,Cp_aortaHTR_smooth*100,'ro-','LineWidth',2)

plot(0:length(trueAIF)-1,trueAIF*100,'k--','LineWidth',2)
hold off
axis([0 160 -0 480])
xlabel('Time (s)')
ylabel('% Signal Change')
legend('True AIF','VS-AIF','HTR-AIF')
paper_setFigureProps;


%%
gfrInd = 8;
dd = 9;
figure
trueKidney = trueCkidney(:,gfrInd);
fitInfo = fitInfoMat{dd,gfrInd};
plot(0:length(trueKidney)-1,trueKidney*100,'k--','LineWidth',2)
hold on
plot(fitInfo.t_aortaVS-1,fitInfo.Ckidney*100,'bo-','LineWidth',2)
hold off

xlabel('Time (s)')
ylabel('% Signal Change')
legend('True Cortex','VS-Cortex')
paper_setFigureProps;

%%
figure
plot(0:length(trueAIF)-1,trueAIF*100,'r-','LineWidth',2)
hold on
plot(0:length(trueAIF)-1,trueCkidney*100,'k--','LineWidth',2)
hold off
axis([0 160 0 450])
xlabel('Time (s)')
ylabel('% Signal Change')
legend('True AIF','True Cortex')
paper_setFigureProps;



%%
% vbVecHTR = zeros(96,1);
% vbVecVS = zeros(96,1);
% cnt = 0;
% for gg = 1:8
%     for dd = 1:12
%         fitInfo = fitInfoMat{dd,gg};
%         cnt = cnt+1;
%         vbVecHTR(cnt) = fitInfo.xFitHTR(3);
%         vbVecVS(cnt) = fitInfo.xFitVS(3);
%     end
% end
% figure
% plot(vbVecHTR,'r-')
% hold on
% plot(vbVecVS,'b-')

%%
% ktransVec = zeros(12,1);
% for gg = 8:8
%     for dd = 1:12
%         fitInfo = fitInfoMat{dd,gg};
%         Ckidney = fitInfo.Ckidney;
%         t_kidney = fitInfo.t_kidney;
%         [xFitTRUE,resnormTRUE] = FitThreeCompartmentAsymmetric(Ckidney,t_kidney,trueAIF,0:length(trueAIF)-1);
%         ktransVec(dd) = xFitTRUE(1);
%     end
% end
% figure
% plot(ktransVec)
%%
fitInfo = fitInfoMat{5,8};
Ckidney = fitInfo.Ckidney;
t_kidney = fitInfo.t_kidney;
trueKidney = trueCkidney(:,8);
t_true = (0:length(trueAIF)-1).';
figure
plot(t_true,trueKidney,'k-')
hold on
plot(t_kidney,Ckidney,'b-')

xFitHTR = fitInfo.xFitHTR;
[xFitTRUE,resnormTRUE] = FitThreeCompartmentAsymmetric(Ckidney,t_kidney,trueAIF,0:length(trueAIF)-1);
xFitIDEAL = [1.3785, 1.50, 0.35, 0.20, 1.8, 0];

xdata = [trueAIF t_true];
[CFitTRUE,Cart,Ctub] = ThreeCompartment(xFitTRUE,xdata);
[CFitIDEAL,Cart,Ctub] = ThreeCompartment(xFitIDEAL,xdata);
[CFitHTR,Cart,Ctub] = ThreeCompartment(xFitHTR,xdata);

figure
plot(t_true,trueKidney,'k-')
hold on
plot(t_true,CFitTRUE,'g-')
% plot(t_true,CFitIDEAL,'r--')
plot(t_true,CFitHTR,'r-')
plot(t_kidney,Ckidney,'b.')



