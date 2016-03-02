
CA

datapath = 'C:\LabviewRoot\Data_files\';
outputpath = 'C:\Users\Wilson Lab\weekly_meetings\04_29_2013_meeting\';
files = {'data_14_07_23_163516.o'};
label = 'brp';
range = {[]};
freqs = [1];
brist = 1;%% if brist = 1, plot bristle
patch = 0; %% if patch = 1, plot patch
output_figure = 0;

ff = 1;
filename = files{ff};
fs = 10000;
data = ReadAllDataFromFile([datapath filename]);
patch_1 = data(2,:)./3; %%patch primary signal
patch_2 = data(4,:)./3*100; %%patch current signal
bristle_1 = data(1,:)./3; %%bristle primary signal
bristle_1 = bristle_1-mean(bristle_1)-27;
bristle_2 = data(3,:)./3*100; %%bristle current signal (1 mv/ pA)
bristle_2 = (bristle_2-mean(bristle_2))/8-50;
% patch_3 = data(6,:)./3*100; %%bristle command potential

scal_const = 10;offset = 55;
piezo_1 = -data(5,:)./3./scal_const;
piezo_1 = (piezo_1-mean(piezo_1));
piezo_2 = data(6,:)./3./scal_const;
piezo_2 = (piezo_2-mean(piezo_2));
patch_1_range = max(patch_1) - min(patch_1);

if ~isempty(range{ff}); patch_1 = patch_1(range{ff});end
if ~isempty(range{ff}); patch_2 = patch_2(range{ff});end
if ~isempty(range{ff}); piezo_1 = piezo_1(range{ff}); piezo_2 = piezo_2(range{ff});end
if ~isempty(range{ff}); bristle_1 = bristle_1(range{ff});end
if ~isempty(range{ff}); bristle_2 = bristle_2(range{ff});end

time_vec = (1:length(patch_1))/fs;
secs = length(patch_1)/fs;
plot_range = 1*fs:secs*fs;

fig1 = figure(1);clf;hold all;set(fig1,'color', 'w', 'position', [1 1 600 800]);
subplot(2,1,1); hold all;box off;
if patch ==1;plot(time_vec, patch_1, 'k', 'linewidth', 0.75);end
if brist ==1; plot(time_vec, bristle_1, 'b', 'linewidth', 0.75);end
plot(time_vec, piezo_1-offset, 'r', 'linewidth', 0.75);
plot(time_vec, piezo_2-offset-5, 'b', 'linewidth', 0.75);

Fs = 10000;x = bristle_1;
N = length(x);
xdft = fft(x);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)).*abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(x):Fs/2;

subplot(2,1,2);
plot(freq,10*log10(psdx)); grid on;
title('Periodogram Using FFT');
xlabel('Frequency (Hz)'); ylabel('Power/Frequency (dB/Hz)');
xlim([0 100]);

% figure;
% [psdestx,Fxx] = periodogram(x,rectwin(length(x)),length(x),Fs);
% plot(Fxx,10*log10(psdestx)); grid on;
% xlabel('Hz'); ylabel('Power/Frequency (dB/Hz)');
% title('Periodogram Power Spectral Density Estimate');
% xlim([0 120]);
