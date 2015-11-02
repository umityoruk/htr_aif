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

alphaVec = 0:0.1:0.6;
A = length(alphaVec);

KtransMatHTR = zeros(D,G);
alphaBest = zeros(D,G);
for gg = 1:G
    for dd = 1:D
        fprintf(['(' num2str(gg) '/' num2str(G) ')' ...
            ' (' num2str(dd) '/' num2str(D) ')\n']);
        
        fitInfo = fitInfoMat{dd,gg};
        
        Ckidney = fitInfo.Ckidney;
        t_kidney = fitInfo.t_kidney;
        
        Cp_aortaHTR = fitInfo.Cp_aortaHTR;
        t_aortaHTR = fitInfo.t_aortaHTR;
                
        
        for aa = 1:A
            alpha = alphaVec(aa);
            if (alpha > 0)
                Cp_aorta_smooth = smooth(Cp_aortaHTR,alpha,'loess').';
            else
                Cp_aorta_smooth = Cp_aortaHTR;
            end
            % Get tail section
            tolerance = 0.05; % Use 0.075 to get the exact figure
            selectLast = 12;
            bb = Cp_aorta_smooth(end-(selectLast-1):end);
            bb = bb(:);
            t_tail = t_aortaHTR(end-(selectLast-1):end);
            AA = [t_tail(:) ones(length(t_tail),1)];
            xx = AA\bb;
            tail = AA*xx;
            tail_max = tail*(1+tolerance);
            tail_min = tail*(1-tolerance);
            if ((sum(bb>tail_max) + sum(bb<tail_min)) == 0)
                % Success
                break;
            end
        end
        [xFitHTR,resnorm] = FitThreeCompartmentAsymmetric(Ckidney,t_kidney,Cp_aorta_smooth,t_aortaHTR);
        alphaBest(dd,gg) = alpha;
        KtransMatHTR(dd,gg) = xFitHTR(1);
        fprintf(['alphaBest = ' num2str(alphaBest(dd,gg)) '\n'])
    end
end
load('targetROI.mat')
GFRMatHTR = KtransMatHTR*targetROI.Vvox*targetROI.RCVoxCnt;
GFRErrHTR = (GFRMatHTR - repmat(GFR,[D 1]))./repmat(GFR,[D 1])*100;
medErrHTR = median(abs(GFRErrHTR),1);
medErrVS = median(abs(GFRErrVS),1);


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




