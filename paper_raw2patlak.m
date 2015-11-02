%% SETUP
paper_loadParameters;
% GFR

load('simResultsRaw.mat')
load('truePSC.mat')
load('targetROI.mat')


fitInfoMat = simResults.fitInfoMat;

trueAIF = truePSC.trueAIF;
trueCkidney = truePSC.Ckidney;
%% APPLY SMOOTHING
G = length(GFR);
D = 12;
KtransMatVS = zeros(D,G);
KtransMatHTR = zeros(D,G);
KtransMatTRUE = zeros(D,G);
KtransMatTRUECX = zeros(D,G);
KtransMatHTRCX = zeros(D,G);
resnormMatVS = zeros(D,G);
resnormMatHTR = zeros(D,G);
resnormMatTRUE = zeros(D,G);
resnormMatTRUECX = zeros(D,G);
resnormMatHTRCX = zeros(D,G);


patlakInfoMat = cell(D,G);
patlakNumMat_VS = zeros(D,G);
patlakNumMat_HTR = zeros(D,G);
for gg = 1:length(GFR)
    for dd = 1:size(fitInfoMat,1)
        fitInfo = fitInfoMat{dd,gg};
        Ckidney = fitInfo.Ckidney;
        t_kidney = fitInfo.t_kidney;
        Cp_aortaVS = fitInfo.Cp_aortaVS;
        t_aortaVS = fitInfo.t_aortaVS;
        Cp_aortaHTR = smooth(fitInfo.Cp_aortaHTR,0.2,'loess');
        CkidneyHTR = smooth(fitInfo.CkidneyHTR,0.2,'loess');
        t_aortaHTR = fitInfo.t_aortaHTR;
        %         [xFitVS,resnormVS] = FitThreeCompartmentAsymmetric(Ckidney,t_kidney,Cp_aortaVS,t_aortaVS);
        %         [xFitHTR,resnormHTR] = FitThreeCompartmentAsymmetric(Ckidney,t_kidney,Cp_aortaHTR,t_aortaHTR);
        %         [xFitTRUE,resnormTRUE] = FitThreeCompartmentAsymmetric(Ckidney,t_kidney,trueAIF,0:length(trueAIF)-1);
        %         [xFitTRUECX,resnormTRUECX] = FitThreeCompartmentAsymmetric(CkidneyHTR,t_aortaHTR,trueAIF,0:length(trueAIF)-1);
        %         [xFitHTRCX,resnormHTRCX] = FitThreeCompartmentAsymmetric(CkidneyHTR,t_aortaHTR,Cp_aortaHTR,t_aortaHTR);
        
        
        
        %% PATLAK
        TRes = 1;
        t_true = 1:TRes:length(trueAIF);
        nPatlakPts=length(Ckidney);
        Cp_aortaVS = interp1(t_aortaVS,Cp_aortaVS,t_true,'linear','extrap');
        Cp_aortaHTR = interp1(t_aortaHTR,Cp_aortaHTR,t_true,'linear','extrap');
        t_kidney = t_aortaVS;
        
        Cp_aortaVS(Cp_aortaVS<0) = 0;
        Cp_aortaHTR(Cp_aortaHTR<0) = 0;
        
        patlakX = zeros(nPatlakPts,1);
        patlakY = zeros(nPatlakPts,1);
        for i=1:nPatlakPts
            % find index of nearest t_true
            [~,indTrue] = min(abs(t_true - t_kidney(i)));
            if (Cp_aortaVS(indTrue)<0.01)
                continue;
            end
            patlakX(i) = simpsonSum(Cp_aortaVS,TRes,1,indTrue)/Cp_aortaVS(indTrue);
            patlakY(i) = Ckidney(i)/Cp_aortaVS(indTrue);
        end
        % sort Patlak points based on patlakX
        patlakMat_VS = sortrows([patlakX patlakY]);
        
        patlakX = zeros(nPatlakPts,1);
        patlakY = zeros(nPatlakPts,1);
        for i=1:nPatlakPts
            % find index of nearest t_true
            [~,indTrue] = min(abs(t_true - t_aortaVS(i)));
            if (Cp_aortaHTR(indTrue)<0.01)
                continue;
            end
            patlakX(i) = simpsonSum(Cp_aortaHTR,TRes,1,indTrue)/Cp_aortaHTR(indTrue);
            patlakY(i) = Ckidney(i)/Cp_aortaHTR(indTrue);
        end
        % sort Patlak points based on patlakX
        patlakMat_HTR = sortrows([patlakX patlakY]);
        
        patlakInfo.patlakMat_VS = patlakMat_VS;
        patlakInfo.patlakMat_HTR = patlakMat_HTR;
        patlakInfoMat{dd,gg} = patlakInfo;
        
        %         figure(100)
        %         plot(patlakMat_HTR(:,1),patlakMat_HTR(:,2),'r.')
        %         hold on
        %         plot(patlakMat_VS(:,1),patlakMat_VS(:,2),'b.')
        %         hold off
        %         pause;
        
        %% BEST FIT TO PATLAK
        ptrangeMin = 10;
        ptrangeMax = 60;
        
        patlakX = patlakMat_VS(:,1);
        patlakY = patlakMat_VS(:,2);
        rangeSelect = find((patlakX >= ptrangeMin) .* (patlakX <= ptrangeMax));
        
        Y = patlakY(rangeSelect);
        X = [patlakX(rangeSelect) ones(length(rangeSelect),1)];
        ab = X\Y;
        patlakNum = ab(1)*60; % [1/min]
        patlakNumMat_VS(dd,gg) = ab_correction(1)*patlakNum + ab_correction(2);
        
%         figure(101)
%         plot(patlakX,patlakY(:),'k*')
%         hold on
%         plot(X(:,1),ab(1)*X(:,1)+ab(2),'r-','LineWidth',2)
%         hold off
%         pause
        
        
        patlakX = patlakMat_HTR(:,1);
        patlakY = patlakMat_HTR(:,2);
        rangeSelect = find((patlakX >= ptrangeMin) .* (patlakX <= ptrangeMax));
        
        Y = patlakY(rangeSelect);
        X = [patlakX(rangeSelect) ones(length(rangeSelect),1)];
        ab = X\Y;
        patlakNum = ab(1)*60; % [1/min]
        patlakNumMat_HTR(dd,gg) = ab_correction(1)*patlakNum + ab_correction(2);
        
%         figure(101)
%         plot(patlakX,patlakY(:),'k*')
%         hold on
%         plot(X(:,1),ab(1)*X(:,1)+ab(2),'r-','LineWidth',2)
%         hold off
%         pause
        
        
        
        %         KtransMatVS(dd,gg) = xFitVS(1);
        %         resnormMatVS(dd,gg) = resnormVS;
        %         KtransMatHTR(dd,gg) = xFitHTR(1);
        %         resnormMatHTR(dd,gg) = resnormHTR;
        %         KtransMatTRUE(dd,gg) = xFitTRUE(1);
        %         resnormMatTRUE(dd,gg) = resnormTRUE;
        %         KtransMatTRUECX(dd,gg) = xFitTRUECX(1);
        %         resnormMatTRUECX(dd,gg) = resnormTRUECX;
        %         KtransMatHTRCX(dd,gg) = xFitHTRCX(1);
        %         resnormMatHTRCX(dd,gg) = resnormHTRCX;
        %
        %         fitInfo.xFitVS = xFitVS;
        %         fitInfo.xFitHTR = xFitHTR;
        %         fitInfo.xFitTRUE = xFitTRUE;
        %         fitInfo.xFitTRUECX = xFitTRUECX;
        %         fitInfo.xFitHTRCX = xFitHTRCX;
        %         fitInfo.Cp_aortaHTR = Cp_aortaHTR;
        %         fitInfo.CkidneyHTR = CkidneyHTR;
        %         fitInfoMat{dd,gg} = fitInfo;
    end
end



GFRMatVS = patlakNumMat_VS*targetROI.Vvox*targetROI.RCVoxCnt;
GFRMatHTR = patlakNumMat_HTR*targetROI.Vvox*targetROI.RCVoxCnt;
% GFRMatTRUE = KtransMatTRUE*targetROI.Vvox*targetROI.RCVoxCnt;
% GFRMatTRUECX = KtransMatTRUECX*targetROI.Vvox*targetROI.RCVoxCnt;
% GFRMatHTRCX = KtransMatHTRCX*targetROI.Vvox*targetROI.RCVoxCnt;
GFRErrVS = (GFRMatVS - repmat(GFR,[D 1]))./repmat(GFR,[D 1])*100;
GFRErrHTR = (GFRMatHTR - repmat(GFR,[D 1]))./repmat(GFR,[D 1])*100;
% GFRErrTRUE = (GFRMatTRUE - repmat(GFR,[D 1]))./repmat(GFR,[D 1])*100;
% GFRErrTRUECX = (GFRMatTRUECX - repmat(GFR,[D 1]))./repmat(GFR,[D 1])*100;
% GFRErrHTRCX = (GFRMatHTRCX - repmat(GFR,[D 1]))./repmat(GFR,[D 1])*100;

medErrVS = median(abs(GFRErrVS),1);
medErrHTR = median(abs(GFRErrHTR),1);

figure
plot(GFR(1:end),medErrVS(1:end),'b-')
hold on
plot(GFR(1:end),medErrHTR(1:end),'r-')
hold off
title('Median Error')
xlabel('GFR [ml/min]')
ylabel('% Error')
legend('VS','HTR')
