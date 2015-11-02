
x = 0:0.5:100;
y = normpdf(x,33,9);

figure
plot(x,y,'k-','LineWidth',2)
hold on
plot([x(x==27) x(x==39)],[y(x==27) y(x==39)],'b*-','LineWidth',2)
plot([x(x==31.5) x(x==34.5)],[y(x==31.5) y(x==34.5)],'r*-','LineWidth',2)
xlabel('Time (s)')
ylabel('Signal Change (\DeltaS)')

%%
x = 0:0.5:100;
y = normpdf(x,10,3);

% figure
plot(x,y,'k-','LineWidth',2)
hold on
plot([x(x==4) x(x==16)],[y(x==4) y(x==16)],'b*-','LineWidth',2)
plot([x(x==8.5) x(x==11.5)],[y(x==8.5) y(x==11.5)],'r*-','LineWidth',2)