CA

% datapath = 'C:\LabviewRoot\Data_files\';
datapath = 'C:\Labviewroot\Data_files\';
% datapath = 'H:\Ephys Data Backup\Datafiles\';
outputpath = 'C:\Users\Wilson Lab\Dropbox\2014_02_25_2014_meeting\';
files = {'data_14_09_22_103212.o'};
label = 'chagal4_uasreachr_attp5_mech';
file_name = files{1};
range = {[]};
freqs = [1];
brist = 1;%% if brist = 1, plot bristle
patch = 0; %% if patch = 1, plot patch
output_figure = 0; %% produce pdf figure
other_piezo =0; %% if 2nd piezo is the primary one, == 1

ff = 1;
filename = files{ff};
fs = 10000;
data = ReadAllDataFromFile([datapath filename]);
patch_1 = data(2,:)./3; %%patch primary signal
patch_2 = data(4,:)./3*100; %%patch current signal
bristle_1 = data(1,:)./3; %%bristle primary signal
bristle_1 = bristle_1-mean(bristle_1);
bristle_2 = data(3,:)./3; %%bristle current signal (1 mv/ pA)
bristle_2 = (bristle_2-mean(bristle_2));
% patch_3 = data(6,:)./3*100; %%bristle command potential


piezo_1 = -data(5,:)./4;
piezo_1 = (piezo_1-mean(piezo_1));
piezo_2 = data(6,:)./4;
piezo_2 = (piezo_2-mean(piezo_2));
% piezo_2 = zeros(1,length(piezo_2));
patch_1_range = max(patch_1) - min(patch_1);

cutoff = 200;%%cutoff frequencies for high-pass filtering data
filts = cutoff/(fs/2);
[x,y] = butter(2,filts,'high');%%bandpass filter between 50 and 200 Hz
filtered_bristle_1 = filter(x, y, bristle_1)';

if other_piezo == 1; piezo = piezo_2; else piezo = piezo_1; end

if ~isempty(range{ff}); patch_1 = patch_1(range{ff});end
if ~isempty(range{ff}); patch_2 = patch_2(range{ff});end
if ~isempty(range{ff}); piezo = piezo(range{ff});end
if ~isempty(range{ff}); piezo_1 = piezo_1(range{ff});end
if ~isempty(range{ff}); piezo_2 = piezo_2(range{ff});end

if ~isempty(range{ff}); bristle_1 = bristle_1(range{ff});end
if ~isempty(range{ff}); bristle_2 = bristle_2(range{ff});end

time_vec = (1:length(patch_1))/fs;
secs = length(patch_1)/fs;
plot_range = 1*fs:secs*fs;

txtfile = [datapath 'pos' filename(5:end-2) '.txt'];
txt_filename = txtfile;
try
values = read_stepper_text_12_2013(txtfile);
x_pos = values(2,:)-mean(values(2,:));x_pos = x_pos-mean(x_pos);
y_pos = values(3,:)-mean(values(3,:));y_pos = y_pos-mean(y_pos);
z_pos = values(4,:)-mean(values(4,:)); z_pos = z_pos-mean(z_pos);
x_pos2 = values(5,:)-mean(values(5,:));x_pos2 = x_pos2-mean(x_pos);
y_pos2 = values(6,:)-mean(values(6,:));y_pos2 = y_pos2-mean(y_pos);
z_pos2 = values(7,:)-mean(values(7,:)); z_pos2 = z_pos2-mean(z_pos);
stepper_failed = [];
catch stepper_failed;
end

stepper_scaling = 10;
offset = -max(patch_1)+2;
piezo_offset = -min(patch_1)+2;
piezo_scaling = (max(piezo)-min(piezo))/3;

fig1 = figure(1);clf;hold all;set(fig1,'color', 'w', 'position', [100 100 800 1000]);
subplot(3,1,1); box off; hold all;
p(1,:) = plot(time_vec, piezo_1/piezo_scaling-piezo_offset, 'r', 'linewidth', 1);
p(2,:) = plot(time_vec, piezo_2/piezo_scaling-piezo_offset, 'b', 'linewidth', 1);
% p(3,:) = plot(time_vec, cond_sig-offset, 'c', 'linewidth', 1);
if isempty(stepper_failed)
p(4,:) = plot(time_vec(1:1000:end),y_pos(1:ceil(length(time_vec)/1000))*stepper_scaling-offset, 'm');
end
if patch ==1;plot(time_vec, patch_1, 'k', 'linewidth', 1);end
if brist ==1; plot(time_vec, bristle_1, 'linewidth', 1, 'color', [0.2 0.2 0.2]);end

% legend([p(1,1) p(2,1) p(3,1) p(4,1)], 'piezo 1', 'piezo 2', 'cond sig', 'stepper', 'location', 'southeast');
% legend('boxoff');
% ylim([-70 0]);
ylabel('mV'); xlabel('time (s)');
title(['entire recording' files], 'interpreter', 'none');

% cutoff = 200;%%cutoff frequencies for high-pass filtering data
% filts = cutoff/(fs/2);
% [x,y] = butter(2,filts,'high');%%bandpass filter between 50 and 200 Hz
% filtered_bristle_1 = filter(x, y, bristle_1)';
% % figure;hold all; plot(bristle_1,'k');plot(filtered_bristle_1,'b');

decs = find(diff(piezo) < -0.25);
incs = find(diff(piezo) > 0.25);
if decs(1) < incs(1);decs(1) = [];end
if any(diff(incs) < 100); incs(diff(incs) < 10) = [];end %%make sure incs are not too close
if any(diff(decs) < 100); decs(diff(decs) < 10) = [];end
if incs(length(incs)) > decs(length(decs)); incs(length(incs)) = [];end

per = 2.*abs(decs(1)-incs(1));
extra = 5000;
patches_1 = [];stims = [];ii = 1;
for ii = 1:length(decs)
if incs(ii)+per < length(patch_1) && incs(ii)-extra > 0
    patches_1(ii,:) = patch_1(incs(ii)-extra:incs(ii)+per)';
    patches_2(ii,:) = patch_2(incs(ii)-extra:incs(ii)+per)';
    stims(ii,:) = piezo(incs(ii)-extra:incs(ii)+per)';
    bristles_1(ii,:) = bristle_1(incs(ii)-extra:incs(ii)+per)';
    bristles_2(ii,:) = bristle_2(incs(ii)-extra:incs(ii)+per)';
ii = ii+1;
end
end


cmap = colormap(autumn(length(decs)));
subplot(3,1,2); hold all;box off;

if patch == 1
for jj = 1:length(decs)
ext = time_vec(1:decs(jj)-incs(jj)+401);
if patch == 1; plot(ext, patch_1((incs(jj)-200):(decs(jj)+200)), 'linewidth', 0.75, 'color', cmap(jj,:));end
plot(ext, piezo((incs(jj)-200):(decs(jj)+200))/(piezo_scaling*2)-piezo_offset, 'b', 'linewidth', 0.75);
end
% ylim([-60 -33]);xlim([0 0.1])
end

if brist == 1
for jj = 1:length(decs)
ext = time_vec(1:decs(jj)-incs(jj)+401);
 plot(ext, filtered_bristle_1((incs(jj)-200):(decs(jj)+200)), 'linewidth', 0.75, 'color', cmap(jj,:))
plot(ext, piezo((incs(jj)-200):(decs(jj)+200))/(piezo_scaling*2)-piezo_offset, 'b', 'linewidth', 0.75);
end
% ylim([-60 -33]);xlim([0 0.1])
end


% subplot(2,2,2); hold all;box off;
% ext = time_vec((incs(2)-200):(decs(2)+200));
% if patch == 1; plot(ext, patch_1((incs(2)-200):(decs(2)+200)), 'k', 'linewidth', 0.75);end
% plot(ext, piezo((incs(2)-200):(decs(2)+200)), 'r', 'linewidth', 0.75);
% if brist ==1; plot(ext, bristle_1((incs(2)-200):(decs(2)+200)), 'b', 'linewidth', 0.75);end;
% % ylim([-70 0]);
% ylabel('mV'); xlabel('time (s)');
% title('1 s of data');
% xlim([time_vec(incs(2)-200) time_vec(decs(2)+200)]);


subplot(3,1,3);hold all;box off;
% confplot((1:size(patches_1,2))/fs,mean(patches_1), std(patches_1), std(patches_1), 'k');hold on;
if brist ==1; confplot((1:size(bristles_1,2))/fs,mean(bristles_1), std(bristles_1), std(bristles_1), 'b');hold on;end
if patch ==1; confplot((1:size(patches_1,2))/fs,mean(patches_1), std(patches_1), std(patches_1), 'k');hold on;end

legs(1) = plot((1:size(patches_1,2))/fs,mean(stims)/piezo_scaling-piezo_offset, 'r', 'linewidth', 0.75);hold on;

% legs(2) = plot((1:size(patches_1,2))/fs,mean(bristles_1), 'b', 'linewidth', 0.75);
% legs(3) = plot((1:size(patches_1,2))/fs,mean(patches_1), 'k', 'linewidth', 0.75);
% legend(legs,'piezo signal', 'bristle signal', 'V_m', 'location', 'northeast');

% ylim([-53 -46]);
% xlim([0.499 0.51])
ylabel('V_m (mv) ');xlabel('time (s)');
title({[filename(1:end-2)], 'average piezo cycle'}, 'interpreter', 'none');
% xlim([1.5 2.5]);
% h =  text(0.52,5, );set(h, 'FontWeight', 'bold', 'Fontsize', 16) 

if output_figure == 1
    try export_fig(fig1,[outputpath '/' filename(1:end-2) '_' label '.pdf'], '-pdf','-nocrop', '-r600' , '-painters', '-rgb');
    catch;end;
end
% % 
