% DO NOT USE THIS FOR THE PAPER
function htrAIF = paper_getHTRAIF2(reconraw,t_htr,AMask,cpiVec,aortaDisco)

showFigures = false;

percentSI_aorta = (aortaDisco - aortaDisco(1))/aortaDisco(1);

roiVoxelCount = sum(sum(sum(AMask)));
roiData = zeros(roiVoxelCount,size(reconraw,4));
for tt=1:size(reconraw,4)
    reconraw3d = reconraw(:,:,:,tt);
    roiData(:,tt) = reconraw3d(AMask);
end
clear reconraw3d


roidata0 = roiData(:,cpiVec==0);
t0 = t_htr(cpiVec==0);


if (showFigures)
    figure
    plot(mean(roidata0,1).')
end



roidata1 = roiData(:,logical( (cpiVec==1) + (cpiVec==4) + (cpiVec==7) ));
cpivec1 = cpiVec(logical( (cpiVec==1) + (cpiVec==4) + (cpiVec==7) ));
t1 = t_htr(logical( (cpiVec==1) + (cpiVec==4) + (cpiVec==7) ));
roidata1_a = roiData(:,cpiVec==1);
roidata1_b = roiData(:,cpiVec==4);
roidata1_c = roiData(:,cpiVec==7);

sig1_a = mean(roidata1_a(:,1:end),2);
sig1_b = mean(roidata1_b(:,1:end),2);
sig1_c = mean(roidata1_c(:,1:end),2);

scl1_a = mean([sig1_a sig1_b sig1_c],2)./sig1_a;
scl1_b = mean([sig1_a sig1_b sig1_c],2)./sig1_b;
scl1_c = mean([sig1_a sig1_b sig1_c],2)./sig1_c;

roidata1(:,cpivec1==1) = roidata1(:,cpivec1==1).*repmat(scl1_a,1,size(roidata1_a,2));
roidata1(:,cpivec1==4) = roidata1(:,cpivec1==4).*repmat(scl1_b,1,size(roidata1_b,2));
roidata1(:,cpivec1==7) = roidata1(:,cpivec1==7).*repmat(scl1_c,1,size(roidata1_c,2));

if (showFigures)
    figure
    plot(mean(roidata1,1).')
end



roidata2 = roiData(:,logical( (cpiVec==2) + (cpiVec==5) + (cpiVec==8) ));
cpivec2 = cpiVec(logical( (cpiVec==2) + (cpiVec==5) + (cpiVec==8) ));
t2 = t_htr(logical( (cpiVec==2) + (cpiVec==5) + (cpiVec==8) ));
roidata2_a = roiData(:,cpiVec==2);
roidata2_b = roiData(:,cpiVec==5);
roidata2_c = roiData(:,cpiVec==8);

sig2_a = mean(roidata2_a(:,1:end),2);
sig2_b = mean(roidata2_b(:,1:end),2);
sig2_c = mean(roidata2_c(:,1:end),2);

scl2_a = mean([sig2_a sig2_b sig2_c],2)./sig2_a;
scl2_b = mean([sig2_a sig2_b sig2_c],2)./sig2_b;
scl2_c = mean([sig2_a sig2_b sig2_c],2)./sig2_c;

roidata2(:,cpivec2==2) = roidata2(:,cpivec2==2).*repmat(scl2_a,1,size(roidata2_a,2));
roidata2(:,cpivec2==5) = roidata2(:,cpivec2==5).*repmat(scl2_b,1,size(roidata2_b,2));
roidata2(:,cpivec2==8) = roidata2(:,cpivec2==8).*repmat(scl2_c,1,size(roidata2_c,2));

if (showFigures)
    figure
    plot(mean(roidata2,1).')
end



roidata3 = roiData(:,logical( (cpiVec==3) + (cpiVec==6) + (cpiVec==9) ));
cpivec3 = cpiVec(logical( (cpiVec==3) + (cpiVec==6) + (cpiVec==9) ));
t3 = t_htr(logical( (cpiVec==3) + (cpiVec==6) + (cpiVec==9) ));
roidata3_a = roiData(:,cpiVec==3);
roidata3_b = roiData(:,cpiVec==6);
roidata3_c = roiData(:,cpiVec==9);

sig3_a = mean(roidata3_a(:,1:end),2);
sig3_b = mean(roidata3_b(:,1:end),2);
sig3_c = mean(roidata3_c(:,1:end),2);

scl3_a = mean([sig3_a sig3_b sig3_c],2)./sig3_a;
scl3_b = mean([sig3_a sig3_b sig3_c],2)./sig3_b;
scl3_c = mean([sig3_a sig3_b sig3_c],2)./sig3_c;

roidata3(:,cpivec3==3) = roidata3(:,cpivec3==3).*repmat(scl3_a,1,size(roidata3_a,2));
roidata3(:,cpivec3==6) = roidata3(:,cpivec3==6).*repmat(scl3_b,1,size(roidata3_b,2));
roidata3(:,cpivec3==9) = roidata3(:,cpivec3==9).*repmat(scl3_c,1,size(roidata3_c,2));

if (showFigures)
    figure
    plot(mean(roidata3,1).')
end

sig0 = mean(roidata0,1).';
sig1 = mean(roidata1,1).';
sig2 = mean(roidata2,1).';
sig3 = mean(roidata3,1).';


percentReconEnd = mean(percentSI_aorta(end-2:end));

start0 = mean(sig0(1:3));
final0 = mean(sig0(end-2:end));
offset0 = (final0 - start0)/percentReconEnd - start0;
adjSig0 = sig0 + offset0;
baselineVal0 = mean(adjSig0(1:3));
percent0 = (adjSig0 - baselineVal0)/baselineVal0;


start1 = mean(sig1(1:3));
final1 = mean(sig1(end-2:end));
offset1 = (final1 - start1)/percentReconEnd - start1;
adjSig1 = sig1 + offset1;
baselineVal1 = mean(adjSig1(1:3));
percent1 = (adjSig1 - baselineVal1)/baselineVal1;


start2 = mean(sig2(1:3));
final2 = mean(sig2(end-2:end));
offset2 = (final2 - start2)/percentReconEnd - start2;
adjSig2 = sig2 + offset2;
baselineVal2 = mean(adjSig2(1:3));
percent2 = (adjSig2 - baselineVal2)/baselineVal2;


start3 = mean(sig3(1:3));
final3 = mean(sig3(end-2:end));
offset3 = (final3 - start3)/percentReconEnd - start3;
adjSig3 = sig3 + offset3;
baselineVal3 = mean(adjSig3(1:3));
percent3 = (adjSig3 - baselineVal3)/baselineVal3;


if (showFigures)
    figure
    plot(t0,percentSI_aorta,'-k','LineWidth',2);
    hold on
    plot(t0,percent0,'b-')
    plot(t1,percent1,'r-')
    plot(t2,percent2,'g-')
    plot(t3,percent3,'c-')
    hold off
end

smooth0 = percent0;
smooth1 = percent1;
smooth2 = percent2;
smooth3 = percent3;



% Find the mix ratio using least squares fit to the tail section
tailBegin = 80;
t_ltr = t0;
C_data = percentSI_aorta(t_ltr >= tailBegin);
A_data = interp1(t0,smooth0,t_ltr(t_ltr >= tailBegin),'linear','extrap');
B1_data = interp1(t1,smooth1,t_ltr(t_ltr >= tailBegin),'linear','extrap');
B2_data = interp1(t2,smooth2,t_ltr(t_ltr >= tailBegin),'linear','extrap');
B3_data = interp1(t3,smooth3,t_ltr(t_ltr >= tailBegin),'linear','extrap');
k = [B1_data-A_data B2_data-A_data B3_data-A_data]\(C_data-A_data);

k0 = 1-sum(k);
k1 = k(1); k2 = k(2); k3 = k(3);

% Mix the estimates using the calculated ratio.
intpercent0 = interp1(t0,smooth0,t_htr,'linear','extrap');
intpercent1 = interp1(t1,smooth1,t_htr,'linear','extrap');
intpercent2 = interp1(t2,smooth2,t_htr,'linear','extrap');
intpercent3 = interp1(t3,smooth3,t_htr,'linear','extrap');
sigAB = k0*intpercent0 + k1*intpercent1 + k2*intpercent2 + k3*intpercent3;

if (showFigures)
    figure
    plot(t0,percentSI_aorta,'-k','LineWidth',2);
    hold on
    plot(t_htr,sigAB,'r-')
    hold off
    percentChangeHTR = sigAB;
end

% % Find the mix ratio using least squares fit to the tail section
% tailBegin = 80;
% t_ltr = t0;
% C_data = percentSI_aorta(t_ltr >= tailBegin);
% A_data = interp1(t0,smooth0,t_ltr(t_ltr >= tailBegin),'linear','extrap');
% B_data = interp1(t1,smooth1,t_ltr(t_ltr >= tailBegin),'linear','extrap');
% k = (A_data-B_data)\(C_data-B_data);
%
%
% % Mix the estimates using the calculated ratio.
% intpercent0 = interp1(t0,smooth0,t_htr,'linear','extrap');
% intpercent1 = interp1(t1,smooth1,t_htr,'linear','extrap');
% sigAB = k*intpercent0 + (1-k)*intpercent1;
%
% figure
% plot(t0,percentSI_aorta,'-k','LineWidth',2);
% hold on
% plot(t_htr,sigAB,'r-')
% hold off
% percentChangeHTR = sigAB;




htrAIF = percentChangeHTR;
end