% Ny = 206;
% Nz = 34;
Ny = 25;
Nz = 20;
AReg = 0.1;
numBs = 3;
splitB = 3;
oneShotAB = false;
accy = 2.5;
accz = 1;
simAcc = true;
patternTemplate = generateDiscoPatternVec(Ny,Nz,AReg,numBs,splitB,oneShotAB,accy,accz,simAcc);
patternTemplateHTR = generateDiscoPatternVec(Ny,Nz,AReg,numBs,splitB,oneShotAB,accy,accz,false);


patA = patternTemplate(:,:,1);
patB1 = sum(patternTemplate(:,:,2:4),3);
patB2 = sum(patternTemplate(:,:,5:7),3);
patB3 = sum(patternTemplate(:,:,8:10),3);

patterns = cat(3,patA,patB1,patB2,patB3);

options = {'kd','ro','gv','bs'};
figure
hold on
for ii = 1:size(patterns,3)
    [row,col] = find(patterns(:,:,ii));
    plot(col,row,options{ii},'MarkerSize',8,'MarkerFaceColor',options{ii}(1));
end
axis([0 Nz+1 0 Ny+1])
legend('A','B_1','B_2','B_3')

