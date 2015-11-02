load('simResults_s0p2.mat')
fitInfoMat_paper = simResults.fitInfoMat;

% load('simResults.mat')
load('simResults_review.mat')
fitInfoMat_review = simResults.fitInfoMat;

[D,G] = size(fitInfoMat_paper);

k21Vec = zeros(D,1);
for gg = 2%1:G
    for dd = 1:D
        fitInfo_paper = fitInfoMat_paper{dd,gg};
        fitInfo_review = fitInfoMat_review{dd,gg};
%         figure(100)
%         plot(fitInfo_paper.Cp_aortaHTR,'b-')
%         hold on
% %         plot(fitInfo_review.Cp_aortaHTR,'r-')
%         Cp_aorta_smooth = smooth(fitInfo_review.Cp_aortaHTR,0.2,'loess').';
%         plot(Cp_aorta_smooth,'r-')
%         hold off
%         pause
        fitInfo = fitInfo_review;
        Ckidney = fitInfo.Ckidney;
        t_kidney = fitInfo.t_kidney;
        
        Cp_aortaHTR = fitInfo.Cp_aortaHTR;
        Cp_aorta_smooth = smooth(Cp_aortaHTR,0.2,'loess').';
%         Cp_aorta_smooth = Cp_aortaHTR;
        t_aortaHTR = fitInfo.t_aortaHTR;
        
        [xFitHTR,resnormHTR] = FitThreeCompartmentAsymmetric(Ckidney,t_kidney,Cp_aorta_smooth,t_aortaHTR);
        k21Vec(dd) = xFitHTR(1);
    end
end
