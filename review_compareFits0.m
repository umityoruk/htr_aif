load('subject_5y.mat')
subjects{1} = subject;
load('subject_9y.mat')
subjects{2} = subject;
load('subject_17y.mat')
subjects{3} = subject;


subId = 1;
figure
hold on
plot(subjects{subId}.t, subjects{subId}.cortex{1}, 'k*')
plot(subjects{subId}.t, subjects{subId}.cortex{2}, 'r*')
plot(subjects{subId}.t_htr, subjects{subId}.fit{1}, 'k-')
plot(subjects{subId}.t_htr, subjects{subId}.fit{2}, 'r-')
legend('left', 'right')

subId = 2;
figure
hold on
plot(subjects{subId}.t, subjects{subId}.cortex{1}, 'k*')
plot(subjects{subId}.t, subjects{subId}.cortex{2}, 'r*')
plot(subjects{subId}.t_htr, subjects{subId}.fit{1}, 'k-')
plot(subjects{subId}.t_htr, subjects{subId}.fit{2}, 'r-')
legend('left', 'right')

subId = 3;
figure
hold on
plot(subjects{subId}.t, subjects{subId}.cortex{1}, 'k*')
plot(subjects{subId}.t, subjects{subId}.cortex{2}, 'r*')
plot(subjects{subId}.t_htr, subjects{subId}.fit{1}, 'k-')
plot(subjects{subId}.t_htr, subjects{subId}.fit{2}, 'r-')
legend('left', 'right')

