% load('simResults_s0p2.mat')
load('simResults_review.mat')
load('truePSC.mat')

paper_loadParameters;
% GFR

GFRErrVS = simResults.GFRErrVS;
GFRErrHTR = simResults.GFRErrHTR;

fitInfoMat = simResults.fitInfoMat;
[D,G] = size(fitInfoMat);
% D = 1;

smoothVec = 0:0.05:0.6;
% smoothVec = 0.2;
medErrAlphaMat = zeros(length(smoothVec),G);
metricAlphaMat = zeros(length(smoothVec),G);
for ss = 1:length(smoothVec)
    fprintf(['(' num2str(ss) '/' num2str(length(smoothVec)) ')\n']);
    alpha = smoothVec(ss);
    KtransMatHTR = zeros(D,G);
    metricMatHTR = zeros(D,G);
    for gg = 1:G
        for dd = 1:D
            fitInfo = fitInfoMat{dd,gg};
            
            Ckidney = fitInfo.Ckidney;
            t_kidney = fitInfo.t_kidney;
            
%             Ckidney = truePSC.Ckidney(:,gg);
%             t_kidney = 1:length(Ckidney);
            
            Cp_aortaHTR = fitInfo.Cp_aortaHTR;
            t_aortaHTR = fitInfo.t_aortaHTR;

%             Cp_aortaHTR = truePSC.trueAIF;
%             t_aortaHTR = 1:length(Cp_aortaHTR);

%             Cp_aortaHTR = fitInfo.Cp_aortaVS;
%             t_aortaHTR = fitInfo.t_aortaVS;
            
            if (alpha > 0)
                Cp_aorta_smooth = smooth(Cp_aortaHTR,alpha,'loess').';
            else
                Cp_aorta_smooth = Cp_aortaHTR;
            end
            [xFitHTR,resnorm] = FitThreeCompartmentAsymmetric(Ckidney,t_kidney,Cp_aorta_smooth,t_aortaHTR);
            KtransMatHTR(dd,gg) = xFitHTR(1);
            metricMatHTR(dd,gg) = resnorm;
        end
    end
    
%     figure
%     imagesc(metricMatHTR.')
%     set(gca, 'YDir', 'normal')
%     ylabel('GFR')
%     xlabel('delay')
    meanMetricHTR = mean(metricMatHTR,1);
    metricAlphaMat(ss,:) = meanMetricHTR;

    load('targetROI.mat')
    GFRMatHTR = KtransMatHTR*targetROI.Vvox*targetROI.RCVoxCnt;
    GFRErrHTR = (GFRMatHTR - repmat(GFR,[D 1]))./repmat(GFR,[D 1])*100;
    medErrHTR = median(abs(GFRErrHTR),1);
    medErrAlphaMat(ss,:) = medErrHTR;
end
medErrVS = median(abs(GFRErrVS),1);

% medErrAlphaMat = metricAlphaMat;

figure
plot(smoothVec, medErrAlphaMat)

figure
[X,Y] = meshgrid(GFR(2:end),smoothVec);
surf(X,Y,medErrAlphaMat(:,2:end))

[XI,YI] = meshgrid( GFR(2):GFR(end),smoothVec(1):0.01:smoothVec(end));
medErrAlphaInterp = interp2(X,Y,medErrAlphaMat(:,2:end),XI,YI);

surf(XI,YI,medErrAlphaInterp)

figure
imagesc(medErrAlphaInterp.')
set(gca, 'YDir', 'normal')
ylabel('GFR')
xlabel('\alpha_L')





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
fitInfoMat = simResults.fitInfoMat;
[D,G] = size(fitInfoMat);
% G = 1;

trueAIF = truePSC.trueAIF;
rmsErrAlphaMatHTR = zeros(length(smoothVec),D);
rmsErrAlphaMatVS = zeros(length(smoothVec),D);
for ss = 1:length(smoothVec)
    fprintf(['(' num2str(ss) '/' num2str(length(smoothVec)) ')\n']);
    alpha = smoothVec(ss);
    rmsErrMatHTR = zeros(D,G);
    rmsErrMatVS = zeros(D,G);
    for gg = 1:G
        for dd = 1:D
            fitInfo = fitInfoMat{dd,gg};
            Cp_aortaHTR = fitInfo.Cp_aortaHTR;
            Cp_aortaVS = fitInfo.Cp_aortaVS;
            t_aortaHTR = fitInfo.t_aortaHTR;
            t_aortaVS = fitInfo.t_aortaVS;
            trueAIF_HTR = interp1(1:length(trueAIF),trueAIF,t_aortaHTR);
            trueAIF_VS = interp1(1:length(trueAIF),trueAIF,t_aortaVS);
            if (alpha > 0)
                Cp_aorta_smooth = smooth(Cp_aortaHTR,alpha,'loess').';
            else
                Cp_aorta_smooth = Cp_aortaHTR;
            end
            rmseHTR = sqrt(mean((trueAIF_HTR(:) - Cp_aorta_smooth(:)).^2));
            rmseVS = sqrt(mean((trueAIF_VS(:) - Cp_aortaVS(:)).^2));
            rmsErrMatHTR(dd,gg) = rmseHTR;
            rmsErrMatVS(dd,gg) = rmseVS;
        end
    end
    rmsErrAlphaMatHTR(ss,:) = mean(rmsErrMatHTR,2);
    rmsErrAlphaMatVS(ss,:) = mean(rmsErrMatVS,2);
end

figure
hold on
plot(smoothVec, mean(rmsErrAlphaMatVS,2),'b-','LineWidth',2)
plot(smoothVec, mean(rmsErrAlphaMatHTR,2),'r-','LineWidth',2)
xlabel('\alpha_L')
ylabel('RMSE')
legend('VS-AIF','HTR-AIF','Location','SouthEast')

figure
[X,Y] = meshgrid(1:D,smoothVec);
surf(X,Y,rmsErrAlphaMatHTR)

[XI,YI] = meshgrid(1:D,smoothVec(1):0.01:smoothVec(end));
rmsErrAlphaInterp = interp2(X,Y,rmsErrAlphaMatHTR,XI,YI);

surf(XI,YI,rmsErrAlphaInterp)

figure
imagesc(rmsErrAlphaInterp.')
set(gca, 'YDir', 'normal')
ylabel('Delay')
xlabel('alpha')

%%

meanErrVS = mean(abs(GFRErrVS),1);
meanErrHTR = mean(abs(GFRErrHTR),1);

maxErrVS = max(abs(GFRErrVS),[],1);
maxErrHTR = max(abs(GFRErrHTR),[],1);

minErrVS = min(abs(GFRErrVS),[],1);
minErrHTR = min(abs(GFRErrHTR),[],1);

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

%%
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
    CkidneyHTR_smooth = fitInfo.CkidneyHTR;
    plot(fitInfo.t_aortaHTR,CkidneyHTR_smooth,'r-','LineWidth',2)
    
    plot(fitInfo.t_aortaVS-1,fitInfo.Ckidney,'b-','LineWidth',2)
%     hold off
%     pause 
end
title('Ckidney estimation over 12 jitters')
xlabel('Time [s]')
ylabel('% Signal Change')
legend('TrueCkidney','VS')

%%



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

medResVS = median(abs(resnormMatVS),1);
medResHTR = median(abs(resnormMatHTR),1);

figure
plot(GFR(2:end),medResVS(2:end),'b-')
hold on
plot(GFR(2:end),medResHTR(2:end),'r-')
hold off
title('Median Resnorm')
xlabel('GFR [ml/min]')
ylabel('Resnorm')
legend('VS','HTR')


%% Scatter plot

GFRMatVS = simResults.GFRMatVS;
GFRMatHTR = simResults.GFRMatHTR;

medGFRVS = median(abs(GFRMatVS),1);
medGFRHTR = median(abs(GFRMatHTR),1);

% Least squares linear fit
A = [GFR(:) ones(length(GFR),1)];
B_VS = medGFRVS(:);
B_HTR = medGFRHTR(:);
ab_VS = A\B_VS;
ab_HTR = A\B_HTR;
medLinearVS = A*ab_VS;
medLinearHTR = A*ab_HTR;

figure
plot(GFR,GFRMatVS(1,:),'bo','markerSize',10)
hold on
plot(GFR,GFRMatHTR(1,:),'rd','markerSize',10)

plot(GFR,GFR,'k--','LineWidth',2)
plot(GFR,GFRMatVS.','bo','markerSize',10)
plot(GFR,GFRMatHTR.','rd','markerSize',10)
% plot(GFR,medGFRVS,'b-','LineWidth',2)
% plot(GFR,medGFRHTR,'r-','LineWidth',2)
plot(GFR,medLinearVS,'b-','LineWidth',2)
plot(GFR,medLinearHTR,'r-','LineWidth',2)
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