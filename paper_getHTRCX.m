function htrAIF = paper_getHTRCX(reconraw,t_htr,CxMask,cpiVec,cxDisco)

showFigures = false;

percentSI_cx = (cxDisco - cxDisco(1))/cxDisco(1);

roiVoxelCount = sum(sum(sum(CxMask)));
roiData = zeros(roiVoxelCount,size(reconraw,4));
for tt=1:size(reconraw,4)
    reconraw3d = reconraw(:,:,:,tt);
    roiData(:,tt) = reconraw3d(CxMask);
end
clear reconraw3d



% Find the peakTime of the contrast
[~,peakTime] = max(sum(roiData,1));

% Assign the preContrast and peakContrast images
preConRoiData = zeros(roiVoxelCount,max(cpiVec)+1);
peakRoiData = zeros(roiVoxelCount,max(cpiVec)+1);
preConIndexVec = [1 2 3 4 5 6 7 8 9 10];
for cpi=0:max(cpiVec)
    cpiInd = cpi + 1;
    %         preConRoiData(:,cpiInd) = roiData(:,cpiInd);
    preConRoiData(:,cpiInd) = roiData(:,preConIndexVec(cpiInd));
    
    % Find nearest time points to the peak time. When there are multiple
    % frames at equal distance, pick the latest one.
    locationReverse = fliplr(find(cpiVec==cpi));
    [~,Ixreverse] = min(abs(locationReverse-peakTime));
    cpiIndPeak = locationReverse(Ixreverse);
    
    peakRoiData(:,cpiInd) = roiData(:,cpiIndPeak);
end

diffRoiData = peakRoiData - preConRoiData;

%% Find the pixels that change the most
topCutoff = 0.25;
cutoffInd = round(size(diffRoiData,1)*topCutoff);
trackIndex = zeros(cutoffInd,size(diffRoiData,2));
for ii=1:size(diffRoiData,2)
    [~,Ix] = sort(diffRoiData(:,ii),'descend');
    IxTop = Ix(1:cutoffInd);
    trackIndex(:,ii) = IxTop;
end


percentChangeMS = zeros(length(cpiVec),1);  % MS = Multi Slice
for ii=1:length(cpiVec)
    cpi = cpiVec(ii);
    cpiInd = cpi + 1;
    
    currentRoiData = roiData(:,ii);
    
    currentMean = mean(currentRoiData(trackIndex(:,cpiInd)));
    
    baseRoiData = preConRoiData(:,cpiInd);
    baseMean = mean(baseRoiData(trackIndex(:,cpiInd)));
    percentChangeMS(ii) = (currentMean-baseMean)./baseMean;
    
    %         % Alternative Method
    %         currentPercentChange = (currentRoiData(trackIndex(:,cpiInd))-baseRoiData(trackIndex(:,cpiInd)))./baseRoiData(trackIndex(:,cpiInd));
    %         percentChangeMS(ii) = mean(currentPercentChange);
end


percentChange1 = percentChangeMS;


% Apply a fix factor for B regions
percentChange1A = percentChange1(cpiVec==0);
Atimes = t_htr(cpiVec==0);
matchLast=3;
percentChange1Fixed = percentChange1;
scalingFactorB = zeros(max(cpiVec),1);
for ii=1:max(cpiVec)
    percentChange1B = percentChange1(cpiVec==ii);
    lastB = percentChange1B(end-matchLast+1:end);
    Btimes = t_htr(cpiVec==ii);
    lastBtimes = Btimes(end-matchLast+1:end);
    % find the nearest A region for each B time
    AInd = zeros(1,matchLast);
    for matchInd = 1:matchLast
        bt = lastBtimes(matchInd);
        [~,AIx] = min(abs(Atimes-bt));
        AInd(matchInd) = AIx;
    end
    lastA = percentChange1A(AInd);
    scalingFactorB(ii) = mean(lastA./lastB);
end

verbose = false;
if (verbose)
    fprintf(['Scaling Factor: \n']);
    disp(scalingFactorB);
end

% scalingThresh = 2;  % times current mean
% if (sum(scalingFactorB > scalingThresh*mean(scalingFactorB)))
%     scalingFactorB(scalingFactorB>scalingThresh*mean(scalingFactorB))...
%         = mean(scalingFactorB(scalingFactorB<=scalingThresh*mean(scalingFactorB)));
%     if (verbose)
%         fprintf(['Scaling factors corrected! New scaling factor: \n']);
%         disp(scalingFactorB);
%     end
% end


for ii=1:max(cpiVec)
    percentChange1Fixed(cpiVec==ii) = percentChange1(cpiVec==ii)*scalingFactorB(ii);
end


%     % DEBUG
%     % Remove the 3rd ring from results
% %     selectionVec = cpiVec==0|cpiVec==1|cpiVec==4|cpiVec==7|cpiVec==2|cpiVec==5|cpiVec==8|cpiVec==3|cpiVec==6|cpiVec==9;
%     selectionVec = cpiVec==0|cpiVec==1|cpiVec==4|cpiVec==7|cpiVec==2|cpiVec==5|cpiVec==8|cpiVec==3|cpiVec==6|cpiVec==9;
%     selectionVec = selectionVec(1:length(tdisco));
%     percentChange1 = percentChange1(selectionVec);
%     percentChange1Fixed = percentChange1Fixed(selectionVec);
%     tdisco = tdisco(selectionVec);
%     cpiVec = cpiVec(selectionVec);
%     % DEBUG_END


t_ltr = t_htr(cpiVec==0);


scaleForComparison = false;
comparisonScaling = 1;
if (scaleForComparison)
    matchLast = 10;
    timeToStartMatch = t_htr(end-matchLast+1);
    comparisonScaling = mean(percentSI_cx(t_ltr>=timeToStartMatch))/mean(percentChange1Fixed(end-matchLast+1:end));
end


%% Alternative Scaling
% [~,discoTailStartInd] = max(percentSI_cx);
% discoTailStart = t_ltr(discoTailStartInd);
discoTailStart = 1;

discoInterp = interp1(t_ltr,percentSI_cx,t_htr,'cubic','extrap');
discoTail = discoInterp(t_htr >= discoTailStart);
cpiVecTail = cpiVec(t_htr >= discoTailStart);
percentChange1Tail = percentChange1(t_htr >= discoTailStart);
scalingFactorAB = zeros(length(unique(cpiVec)),1);
for ii=0:max(cpiVec)
    Bvec = percentChange1Tail(cpiVecTail==ii);
    Dvec = discoTail(cpiVecTail==ii);
    scalingFactorAB(ii+1) = Bvec\Dvec;
    percentChange1Fixed(cpiVec==ii) = percentChange1(cpiVec==ii)*scalingFactorAB(ii+1);
end
comparisonScaling = 1;
%%


% percentChangeSmooth = smooth(percentChange1Fixed,0.075,'loess');
percentChangeSmooth = smooth(percentChange1Fixed,2/1*0.075,'loess');

if (showFigures)
    figure
    plot(t_htr,percentChange1*comparisonScaling,'r-*')
    hold on
    plot(t_htr(cpiVec==0),percentChange1(cpiVec==0)*comparisonScaling,'b*')
    plot(t_ltr,percentSI_cx,'k-')
    hold off
    title('Uncorrected Echo')
    
    figure
    plot(t_htr,percentChange1Fixed*comparisonScaling,'r-*')
    hold on
    plot(t_htr(cpiVec==0),percentChange1Fixed(cpiVec==0)*comparisonScaling,'b*')
    plot(t_ltr,percentSI_cx,'k-')
    hold off
    title('Corrected Echo')
    
    figure
    plot(t_ltr,percentSI_cx,'k-*')
    hold on
    plot(t_htr,percentChangeSmooth*comparisonScaling,'r-*')
    hold off
    title('Smooth Echo')
    
    figure
    plot(t_ltr,percentSI_cx,'k-','LineWidth',2)
    hold on
    plot(t_htr,percentChangeSmooth*comparisonScaling,'r-','LineWidth',2)
    hold off
    title('Smooth Echo')
end


%%
% percentChangeAvg = mean(percentChangeEcho,1);
% comparisonScaling = 1;
% if (scaleForComparison)
%     timeToStartMatch = tdisco(end-matchLast+1);
%     comparisonScaling = mean(percentSI_aorta(t>=timeToStartMatch))/mean(percentChangeAvg(end-matchLast+1:end));
% end
%
%
% figure
% plot(tdisco,percentChangeAvg*comparisonScaling,'r-*')
% hold on
% plot(tdisco(cpiVec==0),percentChangeAvg(cpiVec==0)*comparisonScaling,'b*')
% plot(t,percentSI_aorta,'k-*')
% hold off
% title('Combined Echo')
%
%
% percentChangeSmooth = smooth(percentChangeAvg,0.1,'loess');
%
% figure
% plot(t,percentSI_aorta,'k-*','LineWidth',1)
% hold on
% plot(tdisco,percentChangeSmooth*comparisonScaling,'r-*','LineWidth',1)
% hold off
% legend('AIF-DISCO','AIF-HTR')
% title(['Smooth Combined Echo'])
%
% figure
% plot(t,percentSI_aorta,'k-','LineWidth',2)
% hold on
% plot(tdisco,percentChangeSmooth*comparisonScaling,'r-','LineWidth',2)
% hold off
% legend('AIF-DISCO','AIF-HTR')
% title(['Smooth Combined Echo'])


percentChangeHTR = percentChangeSmooth*comparisonScaling;
percentChangeHTR = percentChange1Fixed;


htrAIF = percentChangeHTR;

end