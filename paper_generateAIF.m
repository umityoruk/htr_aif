%% Initialize
x = [1.459;2.500;0.30046;0.695;0.1200;0.212;4.650;0.0745;13.078;0.953;20/60];

% Best fit to slow injection
% x = [1.459;2.500;0.30046;0.695;0.1200;0.212;4.650;0.0745;13.078;0.953;20/60];
% Sharp Peak version
% x = [1.959;2.500;0.30046;0.695;0.1000;0.212;4.650;0.0745;13.078;0.953;20/60];

%% Adjust parameters
x_old = x;
x = [1.459;2.500;0.30046;0.695;0.1200;0.212;4.650;0.0745;13.078;0.953;20/60];
% x = [1.959;2.500;0.30046;0.695;0.1000;0.212;4.650;0.0745;13.078;0.953;20/60];
tpar = (0:1:200)';
Cp_old = ParametricAIF(x_old,tpar');
Cp = ParametricAIF(x,tpar');

figure
plot(tpar,Cp,'r-','LineWidth',2);
hold on
plot(tpar,Cp_old,'r--')
legend('NewAIF', 'OldAIF')


%% Save the AIF
saveAIF = false;
if (saveAIF)
    trueAIF = Cp;
    t_AIF = tpar;
    save('trueAIF.mat','trueAIF','t_AIF');
end