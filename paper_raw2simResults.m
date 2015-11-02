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
for gg = 1:length(GFR)
    for dd = 1:size(fitInfoMat,1)
        fitInfo = fitInfoMat{dd,gg};
        Ckidney = fitInfo.Ckidney;
        t_kidney = fitInfo.t_kidney;
        Cp_aortaVS = fitInfo.Cp_aortaVS;
        t_aortaVS = fitInfo.t_aortaVS;
        Cp_aortaHTR = smooth(fitInfo.Cp_aortaHTR,0.15,'loess');
        CkidneyHTR = smooth(fitInfo.CkidneyHTR,0.15,'loess');
        t_aortaHTR = fitInfo.t_aortaHTR;
        [xFitVS,resnormVS] = FitThreeCompartmentAsymmetric(Ckidney,t_kidney,Cp_aortaVS,t_aortaVS);
        [xFitHTR,resnormHTR] = FitThreeCompartmentAsymmetric(Ckidney,t_kidney,Cp_aortaHTR,t_aortaHTR);
        [xFitTRUE,resnormTRUE] = FitThreeCompartmentAsymmetric(Ckidney,t_kidney,trueAIF,0:length(trueAIF)-1);
        [xFitTRUECX,resnormTRUECX] = FitThreeCompartmentAsymmetric(CkidneyHTR,t_aortaHTR,trueAIF,0:length(trueAIF)-1);
        [xFitHTRCX,resnormHTRCX] = FitThreeCompartmentAsymmetric(CkidneyHTR,t_aortaHTR,Cp_aortaHTR,t_aortaHTR);
        
        
        
        KtransMatVS(dd,gg) = xFitVS(1);
        resnormMatVS(dd,gg) = resnormVS;
        KtransMatHTR(dd,gg) = xFitHTR(1);
        resnormMatHTR(dd,gg) = resnormHTR;
        KtransMatTRUE(dd,gg) = xFitTRUE(1);
        resnormMatTRUE(dd,gg) = resnormTRUE;
        KtransMatTRUECX(dd,gg) = xFitTRUECX(1);
        resnormMatTRUECX(dd,gg) = resnormTRUECX;
        KtransMatHTRCX(dd,gg) = xFitHTRCX(1);
        resnormMatHTRCX(dd,gg) = resnormHTRCX;
        
        fitInfo.xFitVS = xFitVS;
        fitInfo.xFitHTR = xFitHTR;
        fitInfo.xFitTRUE = xFitTRUE;
        fitInfo.xFitTRUECX = xFitTRUECX;
        fitInfo.xFitHTRCX = xFitHTRCX;
        fitInfo.Cp_aortaHTR = Cp_aortaHTR;
        fitInfo.CkidneyHTR = CkidneyHTR;
        fitInfoMat{dd,gg} = fitInfo;
    end
end



GFRMatVS = KtransMatVS*targetROI.Vvox*targetROI.RCVoxCnt;
GFRMatHTR = KtransMatHTR*targetROI.Vvox*targetROI.RCVoxCnt;
GFRMatTRUE = KtransMatTRUE*targetROI.Vvox*targetROI.RCVoxCnt;
GFRMatTRUECX = KtransMatTRUECX*targetROI.Vvox*targetROI.RCVoxCnt;
GFRMatHTRCX = KtransMatHTRCX*targetROI.Vvox*targetROI.RCVoxCnt;
GFRErrVS = (GFRMatVS - repmat(GFR,[D 1]))./repmat(GFR,[D 1])*100;
GFRErrHTR = (GFRMatHTR - repmat(GFR,[D 1]))./repmat(GFR,[D 1])*100;
GFRErrTRUE = (GFRMatTRUE - repmat(GFR,[D 1]))./repmat(GFR,[D 1])*100;
GFRErrTRUECX = (GFRMatTRUECX - repmat(GFR,[D 1]))./repmat(GFR,[D 1])*100;
GFRErrHTRCX = (GFRMatHTRCX - repmat(GFR,[D 1]))./repmat(GFR,[D 1])*100;

simResults.KtransMatVS = KtransMatVS;
simResults.KtransMatHTR = KtransMatHTR;
simResults.KtransMatTRUE = KtransMatTRUE;
simResults.KtransMatTRUECX = KtransMatTRUECX;
simResults.KtransMatHTRCX = KtransMatHTRCX;
simResults.resnormMatVS = resnormMatVS;
simResults.resnormMatHTR = resnormMatHTR;
simResults.resnormMatTRUE = resnormMatTRUE;
simResults.resnormMatTRUECX = resnormMatTRUECX;
simResults.resnormMatHTRCX = resnormMatHTRCX;
simResults.GFRMatVS = GFRMatVS;
simResults.GFRMatHTR = GFRMatHTR;
simResults.GFRMatTRUE = GFRMatTRUE;
simResults.GFRMatTRUECX = GFRMatTRUECX;
simResults.GFRMatHTRCX = GFRMatHTRCX;
simResults.GFRErrVS = GFRErrVS;
simResults.GFRErrHTR = GFRErrHTR;
simResults.GFRErrTRUE = GFRErrTRUE;
simResults.GFRErrTRUECX = GFRErrTRUECX;
simResults.GFRErrHTRCX = GFRErrHTRCX;
simResults.fitInfoMat = fitInfoMat;
save('simResults_s0p15.mat','simResults');