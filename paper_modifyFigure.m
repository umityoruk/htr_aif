open('fig_signalChange.fig')

h = gcf;

axesObjs = get(h, 'Children'); %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes

objTypes = get(dataObjs, 'Type');  %type of low-level graphics object


objInd = 1;
xdata = get(dataObjs(objInd), 'XData');  %data from low-level grahics objects
ydata = get(dataObjs(objInd), 'YData');
data_ltr = [xdata(:) ydata(:)];

objInd = 3;
xdata = get(dataObjs(objInd), 'XData');  %data from low-level grahics objects
ydata = get(dataObjs(objInd), 'YData');
data_htr = [xdata(:) ydata(:)];


T = size(data_htr,1);
preCon = 0:9;
postCon = [0 1 2 3 0 4 5 6 0 7 8 9];
numReps = ceil((T - length(preCon))/length(postCon));
cpiVec = [preCon repmat(postCon,1,numReps)];
cpiVec = cpiVec(1:T);

t_htr = data_htr(:,1);
t_ltr = data_ltr(:,1);
y_htr = data_htr(:,2);
y_ltr = data_ltr(:,2);

figure
hold on
ms = 15;
plot(t_htr,y_htr,'-','LineWidth',2,'Color',[0.7 0.1 0.1])
plot(t_ltr,y_ltr,'-','LineWidth',2,'Color',[0.1 0.1 0.7])
plot(t_htr(cpiVec==0),y_htr(cpiVec==0),'ks','LineWidth',1,'MarkerSize',ms,'MarkerFaceColor','k');
plot(t_htr(cpiVec==1),y_htr(cpiVec==1),'rd','LineWidth',1,'MarkerSize',ms,'MarkerFaceColor','r');
plot(t_htr(cpiVec==2),y_htr(cpiVec==2),'ro','LineWidth',1,'MarkerSize',ms,'MarkerFaceColor','r');
plot(t_htr(cpiVec==3),y_htr(cpiVec==3),'rx','LineWidth',2,'MarkerSize',ms,'MarkerFaceColor','r');
plot(t_htr(cpiVec==4),y_htr(cpiVec==4),'gd','LineWidth',1,'MarkerSize',ms,'MarkerFaceColor','g');
plot(t_htr(cpiVec==5),y_htr(cpiVec==5),'go','LineWidth',1,'MarkerSize',ms,'MarkerFaceColor','g');
plot(t_htr(cpiVec==6),y_htr(cpiVec==6),'gx','LineWidth',2,'MarkerSize',ms,'MarkerFaceColor','g');
plot(t_htr(cpiVec==7),y_htr(cpiVec==7),'bd','LineWidth',1,'MarkerSize',ms,'MarkerFaceColor','b');
plot(t_htr(cpiVec==8),y_htr(cpiVec==8),'bo','LineWidth',1,'MarkerSize',ms,'MarkerFaceColor','b');
plot(t_htr(cpiVec==9),y_htr(cpiVec==9),'bx','LineWidth',2,'MarkerSize',ms,'MarkerFaceColor','b');




%%
open('fig_normalized.fig')
h = gcf;

axesObjs = get(h, 'Children'); %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes

objTypes = get(dataObjs, 'Type');  %type of low-level graphics object


objInd = 1;
xdata = get(dataObjs(objInd), 'XData');  %data from low-level grahics objects
ydata = get(dataObjs(objInd), 'YData');
data_ltr = [xdata(:) ydata(:)];

objInd = 3;
xdata = get(dataObjs(objInd), 'XData');  %data from low-level grahics objects
ydata = get(dataObjs(objInd), 'YData');
data_htr = [xdata(:) ydata(:)];

t_htr = data_htr(:,1);
t_ltr = data_ltr(:,1);
y_htr = data_htr(:,2);
y_ltr = data_ltr(:,2);

figure
hold on
ms = 15;
plot(t_htr,y_htr,'-','LineWidth',2,'Color',[0.7 0.1 0.1])
plot(t_ltr,y_ltr,'-','LineWidth',2,'Color',[0.1 0.1 0.7])
plot(t_htr(cpiVec==0),y_htr(cpiVec==0),'ks','LineWidth',1,'MarkerSize',ms,'MarkerFaceColor','k');
plot(t_htr(cpiVec==1),y_htr(cpiVec==1),'rd','LineWidth',1,'MarkerSize',ms,'MarkerFaceColor','r');
plot(t_htr(cpiVec==2),y_htr(cpiVec==2),'ro','LineWidth',1,'MarkerSize',ms,'MarkerFaceColor','r');
plot(t_htr(cpiVec==3),y_htr(cpiVec==3),'rx','LineWidth',2,'MarkerSize',ms,'MarkerFaceColor','r');
plot(t_htr(cpiVec==4),y_htr(cpiVec==4),'gd','LineWidth',1,'MarkerSize',ms,'MarkerFaceColor','g');
plot(t_htr(cpiVec==5),y_htr(cpiVec==5),'go','LineWidth',1,'MarkerSize',ms,'MarkerFaceColor','g');
plot(t_htr(cpiVec==6),y_htr(cpiVec==6),'gx','LineWidth',2,'MarkerSize',ms,'MarkerFaceColor','g');
plot(t_htr(cpiVec==7),y_htr(cpiVec==7),'bd','LineWidth',1,'MarkerSize',ms,'MarkerFaceColor','b');
plot(t_htr(cpiVec==8),y_htr(cpiVec==8),'bo','LineWidth',1,'MarkerSize',ms,'MarkerFaceColor','b');
plot(t_htr(cpiVec==9),y_htr(cpiVec==9),'bx','LineWidth',2,'MarkerSize',ms,'MarkerFaceColor','b');



%%
open('fig_smooth.fig')
h = gcf;

axesObjs = get(h, 'Children'); %axes handles
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes

objTypes = get(dataObjs, 'Type');  %type of low-level graphics object


objInd = 2;
xdata = get(dataObjs(objInd), 'XData');  %data from low-level grahics objects
ydata = get(dataObjs(objInd), 'YData');
data_ltr = [xdata(:) ydata(:)];

objInd = 1;
xdata = get(dataObjs(objInd), 'XData');  %data from low-level grahics objects
ydata = get(dataObjs(objInd), 'YData');
data_htr = [xdata(:) ydata(:)];

t_htr = data_htr(:,1);
t_ltr = data_ltr(:,1);
y_htr = data_htr(:,2);
y_ltr = data_ltr(:,2);

figure
hold on
plot(t_ltr,y_ltr,'-','LineWidth',2,'Color',[0.1 0.1 0.7])
plot(t_htr,y_htr,'-','LineWidth',2,'Color',[0.7 0.1 0.1])
