clear all;

file = '2016_02_05_0015.abf';

ephys_dir = 'E:\ephys_data\';
[d,si,h]=abfload([ephys_dir file]);
fs = si/0.01;
voltage= d(:,1);
current = d(:,2);
piezo_cmd = d(:,3);
piezo_sns = d(:,4);
bristle = d(:,5);

figure; hold all;
% plot((1:length(voltage))/fs,voltage);
plot((1:length(voltage))/fs,bristle);

% plot((1:length(voltage))/fs,piezo_cmd);
plot((1:length(voltage))/fs,piezo_sns-65);
xlabel('time (s)');