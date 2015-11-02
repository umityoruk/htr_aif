topCutoffVec = 0.05:0.1:0.95;
% topCutoffVec = 0.25;
load('simResults_review_25.mat')
load('truePSC.mat')
load('targetROI.mat')

paper_loadParameters;
% GFR

GFRErrVS = simResults.GFRErrVS;
GFRErrHTR = simResults.GFRErrHTR;

fitInfoMat = simResults.fitInfoMat;
[D,G] = size(fitInfoMat);
% D = 1;

medErrAlphaMat = [];
for cc = 1:length(topCutoffVec)
    fprintf(['(' num2str(cc) '/' num2str(length(topCutoffVec)) ')\n']);
    topCutoff = topCutoffVec(cc);
    load(['simResults_review_' num2str(topCutoff*100) '.mat']);
    fitInfoMat = simResults.fitInfoMat;
    [D,G] = size(fitInfoMat);
    KtransMatHTR = zeros(D,G);
    metricMatHTR = zeros(D,G);
    if (isempty(medErrAlphaMat))
        medErrAlphaMat = zeros(length(topCutoffVec),G);
    end
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
            

            Cp_aorta_smooth = smooth(Cp_aortaHTR,0.2,'loess').';

            [xFitHTR,resnorm] = FitThreeCompartmentAsymmetric(Ckidney,t_kidney,Cp_aorta_smooth,t_aortaHTR);
            KtransMatHTR(dd,gg) = xFitHTR(1);
            metricMatHTR(dd,gg) = resnorm;
        end
    end

    GFRMatHTR = KtransMatHTR*targetROI.Vvox*targetROI.RCVoxCnt;
    GFRErrHTR = (GFRMatHTR - repmat(GFR,[D 1]))./repmat(GFR,[D 1])*100;
    medErrHTR = median(abs(GFRErrHTR),1);
    medErrAlphaMat(cc,:) = medErrHTR;
end


save('medErrAlphaMat.mat', 'medErrAlphaMat');
%%
load('medErrAlphaMat.mat');


figure
[X,Y] = meshgrid(GFR(2:end),topCutoffVec);
surf(X,Y,medErrAlphaMat(:,2:end))

[XI,YI] = meshgrid( GFR(2):0.1:GFR(end),topCutoffVec(1):0.01:topCutoffVec(end));
medErrAlphaInterp = interp2(X,Y,medErrAlphaMat(:,2:end),XI,YI);

figure
surf(XI,YI,medErrAlphaInterp)

figure
imagesc(YI(:,1),XI(1,:),medErrAlphaInterp.')
set(gca, 'YDir', 'normal')
ylabel('GFR')
xlabel('\alpha')


%%
medErrHTR = medErrAlphaMat(10,:);

GFRStart = 2;
GFR_line = GFR(GFRStart):GFR(end);
% medErrVS_line = interp1(GFR(GFRStart:end),medErrVS(GFRStart:end),GFR_line,'cubic');
medErrHTR_line = interp1(GFR(GFRStart:end),medErrHTR(GFRStart:end),GFR_line,'cubic');

figure
% plot(GFR(GFRStart:end),medErrVS(GFRStart:end),'bo','markerSize',10,'LineWidth',2,'MarkerFaceColor','w')
hold on
plot(GFR(GFRStart:end),medErrHTR(GFRStart:end),'rd','markerSize',10,'LineWidth',2,'MarkerFaceColor','w')
plot(GFR(GFRStart:end),ones(size(GFR(GFRStart:end))).*5,'k--','LineWidth',2)
plot(GFR(GFRStart:end),ones(size(GFR(GFRStart:end))).*10,'k--','LineWidth',2)
plot(GFR(GFRStart:end),ones(size(GFR(GFRStart:end))).*15,'k--','LineWidth',2)
plot(GFR(GFRStart:end),ones(size(GFR(GFRStart:end))).*20,'k--','LineWidth',2)
plot(GFR(GFRStart:end),ones(size(GFR(GFRStart:end))).*25,'k--','LineWidth',2)
plot(GFR(GFRStart:end),ones(size(GFR(GFRStart:end))).*30,'k--','LineWidth',2)
plot(GFR(GFRStart:end),ones(size(GFR(GFRStart:end))).*35,'k--','LineWidth',2)
% plot(GFR_line,medErrVS_line,'b-','LineWidth',2)
plot(GFR_line,medErrHTR_line,'r-','LineWidth',2)

% plot(GFR(GFRStart:end),medErrVS(GFRStart:end),'bo','markerSize',10,'LineWidth',2,'MarkerFaceColor','w')
plot(GFR(GFRStart:end),medErrHTR(GFRStart:end),'rd','markerSize',10,'LineWidth',2,'MarkerFaceColor','w')

hold off
xlabel('GFR (ml/min)')
ylabel('% Error')
paper_setFigureProps;


