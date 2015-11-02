addpath(genpath('./Source'));

%% Load Global Parameters
paper_loadParameters;
% baseFolderName
% hct

%%
load([baseFolderName '/' 'basePhantom.mat'])
load([baseFolderName '/' 'phantomMask.mat'])
load('targetROI.mat')

applyHTR = true;
%%
pad_left = 25;
pad_right = 25;
Ny = size(basePhantom,2)-(pad_left + pad_right);
Nz = size(basePhantom,3);
AReg = 0.1;
numBs = 3;
splitB = 3;
oneShotAB = false;
accy = 2.5;
accz = 1;
simAcc = true;
patternTemplate = generateDiscoPatternVec(Ny,Nz,AReg,numBs,splitB,oneShotAB,accy,accz,simAcc);
patternTemplateHTR = generateDiscoPatternVec(Ny,Nz,AReg,numBs,splitB,oneShotAB,accy,accz,false);
% simple4DViewer(patternTemplate)

%%
GFR = targetROI.GFR;
% T = size(basePhantom,4);
T = tmax+1;
G = length(GFR);
D = 12;
respRepeat = 3; % 1s acq + 2s pause
KtransMatVS = zeros(D,G);
KtransMatHTR = zeros(D,G);
KtransMatHTRCX = zeros(D,G);
resnormMatVS = zeros(D,G);
resnormMatHTR = zeros(D,G);
resnormMatHTRCX = zeros(D,G);
fitInfoMat = cell(D,G);
for gg = 1:G
    newPhantom = basePhantom(:,:,:,1:T);
    for tt = 1:T
        newPhantom(:,:,:,tt) = paper_buildPhantomAtTime(basePhantom,phantomMask,targetROI,gg,tt);
    end
    
    for dd= 1:D
        fprintf(['(' num2str(gg) '/' num2str(G) ',' ...
            num2str(dd) '/' num2str(D) ')\n']);
        
        cpiVec_pre = 0:9;
        cpiVec_post = ones(1,ceil(T-10/respRepeat)*respRepeat)*-1;
        cpiTemplate = [0 1 2 3 0 4 5 6 0 7 8 9];
        for tt = 1:(length(cpiVec_post)/respRepeat)
            cpiVec_post((tt-1)*3+1) = cpiTemplate(mod(tt-1,length(cpiTemplate))+1);
        end
        cpiVec_post = circshift(cpiVec_post,[0 dd-1]);
        cpiVec = [cpiVec_pre cpiVec_post];
        cpiVec = cpiVec(1:T);
        
        timePointsSampled = find(cpiVec~=-1);
        cpiVecCondensed = cpiVec(cpiVec~=-1);
        
        patternVec = zeros(size(patternTemplate,1),size(patternTemplate,2),length(cpiVecCondensed));
        patternVecHTR = zeros(size(patternTemplateHTR,1),size(patternTemplateHTR,2),length(cpiVecCondensed));
        for tt = 1:length(cpiVecCondensed)
            patternVec(:,:,tt) = patternTemplate(:,:,cpiVecCondensed(tt)+1);
            patternVecHTR(:,:,tt) = patternTemplateHTR(:,:,cpiVecCondensed(tt)+1);
        end
        phantom = newPhantom(:,:,:,timePointsSampled);
        phantomInfo = generatePhantomInfo(phantom,timePointsSampled,pad_left,pad_right);
        samplingInfo.cpiVec = cpiVecCondensed;
        samplingInfo.tCpi = timePointsSampled;
        samplingInfo.parameters.oneShotAB = false;
        viewshareInfo = generateDiscoViewshareInfo(samplingInfo);
        t_htr = timePointsSampled;
        
        ksp = advancedSampleUsingPatternVec(phantom,phantomInfo,patternVec);
        [ksp_vs, t_vs] = advancedViewshare(ksp, viewshareInfo);
        
        recon = advancedRecon(ksp_vs, phantomInfo);
        
        
        aortaNew = averageFromMask(recon,phantomMask.AMask);
        cortexNew = averageFromMask(recon,phantomMask.RCMask);
        
        Ckidney = (cortexNew - cortexNew(1))/cortexNew(1);
        t_kidney = t_vs.';
        
        ksp = advancedSampleUsingPatternVec(phantom,phantomInfo,patternVecHTR);
        recon2 = advancedRecon(ksp, phantomInfo);


        t_aortaHTR = t_htr.';
        Cp_aortaHTR = paper_getHTRAIF(recon2,t_aortaHTR,phantomMask.AMask,cpiVecCondensed,aortaNew)/(1-hct);
        t_aortaVS = t_vs.';
        Cp_aortaVS = ((aortaNew - aortaNew(1))/aortaNew(1))/(1-hct);
        
        CkidneyHTR = paper_getHTRCX(recon2,t_aortaHTR,phantomMask.RCMask,cpiVecCondensed,cortexNew);

        
%         figure
%         plot(0:length(trueAIF)-1,trueAIF,'k-','LineWidth',2);
%         hold on
%         plot(t_aortaHTR,smooth(Cp_aortaHTR,0.2,'loess'),'r-')
%         plot(t_aortaVS,Cp_aortaVS,'b-')
% %         plot(t_aortaHTR,Cp_aortaHTR,'r-')



%         %% DEBUG CODE (DETERMINE SENSITIVITY TO LOESS)
%         smoothVec = 0:0.005:0.4;
%         rmseVec = ones(length(smoothVec),1);
%         for ss = 1:length(smoothVec)
%             alpha = smoothVec(ss);
%             if (alpha > 0)
%                 Cp_aortaHTR_smooth = smooth(Cp_aortaHTR,alpha,'loess');
%             else
%                 Cp_aortaHTR_smooth = Cp_aortaHTR;
%             end
%             trueAIF_htr = interp1(0:length(trueAIF)-1, trueAIF, t_aortaHTR);
%             rmseVec(ss) = sqrt(mean((Cp_aortaHTR_smooth - trueAIF_htr).^2));
%         end
%         
%         figure
%         plot(smoothVec,rmseVec)
        
        %%

% 
%         figure
%         load('truePSC.mat')
%         trueKidney = truePSC.Ckidney(:,gg);
%         plot(t_aortaHTR,CkidneyHTR,'r-')
%         hold on
%         plot(t_aortaVS,Ckidney,'b-')
%         plot(0:length(trueKidney)-1,trueKidney,'k-')
        
        % Three Compartment Fit
        [xFitVS,resnormVS] = FitThreeCompartmentAsymmetric(Ckidney,t_kidney,Cp_aortaVS,t_aortaVS);
        [xFitHTR,resnormHTR] = FitThreeCompartmentAsymmetric(Ckidney,t_kidney,Cp_aortaHTR,t_aortaHTR);
        [xFitHTRCX,resnormHTRCX] = FitThreeCompartmentAsymmetric(CkidneyHTR,t_aortaHTR,Cp_aortaHTR,t_aortaHTR);
        
        %% TEST CODE
        smoothVec = 0.025:0.025:0.2;
        resnormVec = ones(length(smoothVec),1);
        for ss = 1:length(smoothVec)
            Cp_aorta_smooth = smooth(Cp_aortaHTR,smoothVec(ss),'loess').';
            ydata = Ckidney;
            xdata = [Cp_aorta_smooth,t_htr];
            [xFit,resnorm] = FitThreeCompartmentAsymmetric(Ckidney,t_kidney,Cp_aorta_smooth,t_aortaHTR);
            if (sum(resnormVec < resnorm) == 0)
                xFitBest = xFit;
                resnormBest = resnorm;
                smoothBest = smoothVec(ss);
            end
            resnormVec(ss) = resnorm;
        end
        
        xFitHTR = xFitBest;
        resnormHTR = resnormBest;
        fitInfo.smoothBest = smoothBest;
        %%
        
        KtransMatVS(dd,gg) = xFitVS(1);
        resnormMatVS(dd,gg) = resnormVS;
        KtransMatHTR(dd,gg) = xFitHTR(1);
        resnormMatHTR(dd,gg) = resnormHTR;
        KtransMatHTRCX(dd,gg) = xFitHTRCX(1);
        resnormMatHTRCX(dd,gg) = resnormHTRCX;
        
        fitInfo.xFitVS = xFitVS;
        fitInfo.xFitHTR = xFitHTR;
        fitInfo.xFitHTRCX = xFitHTRCX;
        fitInfo.cxSig = cortexNew;
        fitInfo.aortaSig = aortaNew;
        fitInfo.Ckidney = Ckidney;
        fitInfo.CkidneyHTR = CkidneyHTR;
        fitInfo.t_kidney = t_kidney;
        fitInfo.Cp_aortaVS = Cp_aortaVS;
        fitInfo.Cp_aortaHTR = Cp_aortaHTR;
        fitInfo.t_aortaVS = t_aortaVS;
        fitInfo.t_aortaHTR = t_aortaHTR;
        fitInfoMat{dd,gg} = fitInfo;
        
        
        
        
%         % Display the fit
%         xdata = [Cp_aorta t_aorta];
%         xFit2 = xFit;
%         xFit2(1) = 0.1723;
%         xFit2 = [k21Vec_R(1) k12 vb d tau];
%         xdata = [trueAIF (0:T-1).']
%         [xFit,resnorm] = FitThreeCompartmentAsymmetric(Ckidney,t_kidney,trueAIF,(0:T-1).');
%         [Cfit, Cart, Ctub] = ThreeCompartment(xFit2,xdata);
%         figure
%         plot(t_kidney,Ckidney,'k-','LineWidth',2)
%         hold on
%         plot(0:T-1,Cfit,'g-')
    end
end

GFRMatVS = KtransMatVS*targetROI.Vvox*targetROI.RCVoxCnt;
GFRMatHTR = KtransMatHTR*targetROI.Vvox*targetROI.RCVoxCnt;
GFRMatHTRCX = KtransMatHTRCX*targetROI.Vvox*targetROI.RCVoxCnt;
GFRErrVS = (GFRMatVS - repmat(GFR,[D 1]))./repmat(GFR,[D 1])*100;
GFRErrHTR = (GFRMatHTR - repmat(GFR,[D 1]))./repmat(GFR,[D 1])*100;
GFRErrHTRCX = (GFRMatHTRCX - repmat(GFR,[D 1]))./repmat(GFR,[D 1])*100;

simResults.KtransMatVS = KtransMatVS;
simResults.KtransMatHTR = KtransMatHTR;
simResults.KtransMatHTRCX = KtransMatHTRCX;
simResults.resnormMatVS = resnormMatVS;
simResults.resnormMatHTR = resnormMatHTR;
simResults.resnormMatHTRCX = resnormMatHTRCX;
simResults.GFRMatVS = GFRMatVS;
simResults.GFRMatHTR = GFRMatHTR;
simResults.GFRMatHTRCX = GFRMatHTRCX;
simResults.GFRErrVS = GFRErrVS;
simResults.GFRErrHTR = GFRErrHTR;
simResults.GFRErrHTRCX = GFRErrHTRCX;
simResults.fitInfoMat = fitInfoMat;
save('simResults.mat','simResults');