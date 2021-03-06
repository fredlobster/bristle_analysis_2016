CA

datapath = 'C:\LabviewRoot\Data_files\';
spikepath = 'C:\LabviewRoot\Data_files\spikes\';
outputpath = 'C:\Users\Wilson Lab\Dropbox\2014_03_03_meeting\';
files = {'data_14_11_13_161150.o'};
range = {[]};
output_figure = 1;
other_piezo = 0; %% if 2nd piezo is the primary one, == 1
spike_detect = 1;%% 1 if new bristle spikes are to be detected, 0 if not
brist = 1;%% if brist = 1, plot bristle
patch = 0; %% if patch = 1, plot patch

% files = {'data_21_09_12_130320.o'};
% ranges = {[1:3.3e5]};
% freqs = [0.8];

% for ff = 1: length(files)
ff = 1;
filename = files{ff};
fs = 10000;
data = ReadAllDataFromFile([datapath filename]);
spike_file = [spikepath filename(1:end-2) '_spikes.mat'];
patch_1 = data(2,:)./3; %%patch primary signal
patch_2 = data(4,:)./3*100; %%patch current signal
bristle_1 = data(1,:)./3; %%bristle primary signal
bristle_1 = bristle_1-mean(bristle_1)-35;
bristle_2 = data(3,:)./3*100; %%bristle current signal (1 mv/ pA)
bristle_2 = (bristle_2-mean(bristle_2))/8-50;
% patch_3 = data(6,:)./3*100; %%bristle command potential

piezo_1 = -data(5,:)./4;
piezo_1 = (piezo_1-mean(piezo_1));
piezo_2 = data(6,:)./4;
piezo_2 = (piezo_2-mean(piezo_2));
% piezo_2 = zeros(1,length(piezo_2));
patch_1_range = max(patch_1) - min(patch_1);

if ~isempty(range{ff}); patch_1 = patch_1(range{ff});end
if ~isempty(range{ff}); patch_2 = patch_2(range{ff});end
if ~isempty(range{ff}); piezo_1 = piezo_1(range{ff});end
if ~isempty(range{ff}); piezo_2 = piezo_2(range{ff});end
if ~isempty(range{ff}); bristle_1 = bristle_1(range{ff});end
if ~isempty(range{ff}); bristle_2 = bristle_2(range{ff});end

cutoff = 200;%%cutoff frequencies for high-pass filtering data
filts = cutoff/(fs/2);
[x,y] = butter(2,filts,'high');%%bandpass filter between 50 and 200 Hz
filtered_bristle_1 = filter(x, y, bristle_1)';

if other_piezo == 1; piezo = piezo_2; else piezo = piezo_1; end

if exist(spike_file);load(spike_file);spike_detect =0;end

time_vec = (1:length(patch_1))/fs;
secs = length(patch_1)/fs;
plot_range = 1*fs:secs*fs;

%% if necessary, pull out spikes
if spike_detect == 1
fly_start = 1;
baf = 0; %%number of samples to store before and after stimulus
fs = 10000;
spikes.present = 1; %%1 if its a spiking cell, empty if it's not
spikes.approval = 1; %%approval (== 0) to run without asking for the spike distance threshold
spikes.approval_time = 4;
spikes.threshold_factor = 1.4; %%threshold over which to find peaks (# of sds)
spikes.sd = 0.25; %%%made up standard deviation of the data (complicated by steady state deflection)
spikes.threshold = spikes.threshold_factor*spikes.sd;
spikes.spike_count_range = 1:length(bristle_1)+baf; %%range from which to extract spikes
spikes.filts = cutoff/(fs/2);
spikes.spikeDist_threshold = 10;
spikes.spikeTemplateWidth = 55; %%number of samples for the spike template
spikes.spikes2avg = 4;

spiketemplate = [];
fs = 10000;
full_length = length(bristle_1);
num_conditions = 1;

for j = 1:num_conditions  % 36 conditions, create an array of the four important variables for each condition
    Spikes(j).bristle_1 = [];  
    Spikes(j).piezo = []; 
    Spikes(j).spike_locs = [];
    Spikes(j).spike_locs_inds = [];
    Spikes(j).num_spikes = [];
    Spikes(j).spike_ts = [];
end

Index = 1;
Spikes(Index).bristle_1    =   [Spikes(Index).bristle_1; bristle_1'];
Spikes(Index).piezo  =   [Spikes(Index).piezo; piezo'];

[spike_locs, spikesBelowThresh, spiketemplate, spikeDist_threshold] = extracellular_spike_extractor_09_2012(filtered_bristle_1(spikes.spike_count_range)', bristle_1(spikes.spike_count_range)',...
    spikes.threshold,spiketemplate, spikes.spikeTemplateWidth, spikes.spikes2avg, spikes.spikeDist_threshold, spikes.approval, spikes.approval_time);
       Spikes(Index).num_spikes =   [Spikes(Index).num_spikes; length(spike_locs)];%%each row is a list of spike locations       
       Spikes(Index).spike_ts =   [Spikes(Index).spike_ts; spikesBelowThresh'];%%each row is a single spike      
       Spikes(Index).spike_locs =   [Spikes(Index).spike_locs; spike_locs];%%the first column gives the spike location
       trial_num = size(Spikes(Index).bristle_1,1);
       Spikes(Index).spike_locs_inds =   [Spikes(Index).spike_locs_inds; trial_num*ones(length(spike_locs),1)];%%the 2nd column gives the trial in which that spike location resides
       
       spikes.spikeDist_threshold = spikeDist_threshold; %%use the new threshold for the remainder of this cell

save([spikepath filename(1:end-2) '_spikes'], 'Spikes');
end

% decs = find(diff(piezo) < -0.25);
% incs = find(diff(piezo) > 0.25);
% if decs(1) < incs(1);decs(1) = [];end
% if any(diff(incs) < 100); incs(diff(incs) < 10) = [];end %%make sure incs are not too close
% if any(diff(decs) < 100); decs(diff(decs) < 10) = [];end
% if incs(length(incs)) > decs(length(decs)); incs(length(incs)) = [];end
% 
% per = 2.*abs(decs(1)-incs(1));
% extra = 5000;
% patches_1 = [];stims = [];ii = 1;
% for ii = 1:length(decs)
% if incs(ii)+per < length(patch_1) && incs(ii)-extra > 0
%     patches_1(ii,:) = patch_1(incs(ii)-extra:incs(ii)+per)';
%     patches_2(ii,:) = patch_2(incs(ii)-extra:incs(ii)+per)';
%     stims(ii,:) = piezo(incs(ii)-extra:incs(ii)+per)';
%     bristles_1(ii,:) = bristle_1(incs(ii)-extra:incs(ii)+per)';
%     bristles_2(ii,:) = bristle_2(incs(ii)-extra:incs(ii)+per)';
% ii = ii+1;
% end
% end

piezo_offset = -min(patch_1)+2;
piezo_scaling = 3;
add_to_bristle = 15;

fig2 = figure(2);clf;hold all;set(fig2,'color', 'w', 'position', [1 1 800 1050]);
% cmap = colormap(autumn(length(decs)));
% cmap2 = colormap(winter(length(decs)));

box off; hold all;title(filename, 'interpreter', 'none');
p(1,:) = plot(time_vec, piezo_1/piezo_scaling-piezo_offset, 'b', 'linewidth', 1);
p(2,:) = plot(time_vec, piezo_2/(2*piezo_scaling)-piezo_offset-15, 'r', 'linewidth', 1);
% p(3,:) = plot(time_vec, cond_sig-offset, 'c', 'linewidth', 1);
if patch ==1;plot(time_vec, patch_1, 'k', 'linewidth', 1);end
if brist ==1; plot(time_vec, bristle_1+add_to_bristle, 'linewidth', 1, 'color', [0.3 0.3 0.3]);end

bottom = -25;

    plot([Spikes.spike_locs Spikes.spike_locs]/fs,...
        [+bottom 0.5+bottom], 'k-');


% subplot(4,1,2); hold all;box off;title('all trials');
% buffer = 200;
% for jj = 1:length(decs)
% ext = time_vec(1:decs(jj)-incs(jj)+2*buffer+1);
%  if brist == 1; plot(ext,add_to_bristle+bristle_1((incs(jj)-buffer):(decs(jj)+buffer)), 'linewidth', 0.75, 'color', cmap2(jj,:));end
% if patch == 1; plot(ext, patch_1((incs(jj)-200):(decs(jj)+200)), 'linewidth', 0.75, 'color', cmap(jj,:));end
% plot(ext, piezo((incs(jj)-buffer):(decs(jj)+buffer))/piezo_scaling-piezo_offset, 'k', 'linewidth', 0.75);
% end
% 
% bottom = -25+add_to_bristle;
% for jj = 1:length(decs)-1
%     belongs = Spikes.spike_locs > incs(jj) & Spikes.spike_locs < decs(jj)+extra;
% if any(belongs)
%     plot([Spikes.spike_locs(belongs)-decs(jj)+extra+buffer Spikes.spike_locs(belongs)-decs(jj)+extra+buffer]/fs,...
%         [-jj+bottom -jj+0.5+bottom], 'k-');
% end
% end
% 
% axis tight;% ylim([-60 -33]);
% xlim([0 0.55]);
% 
% 
% % subplot(4,1,3:4);hold all;box off;title('averages');
% % if brist ==1; confplot((1:size(bristles_1,2))/fs,mean(bristles_1)-10, std(bristles_1), std(bristles_1), 'b');hold on;end
% % if patch ==1; confplot((1:size(patches_1,2))/fs,mean(patches_1), std(patches_1), std(patches_1), 'r');hold on;end
% % plot((1:size(stims,2))/fs,mean(stims)/piezo_scaling-piezo_offset, 'k');
% 
% subplot(4,1,3:4);hold all;box off;title('averages');
% if brist ==1; plot((1:size(bristles_1,2))/fs,mean(bristles_1)-10, 'b');hold on;end
% if patch ==1; plot((1:size(patches_1,2))/fs,mean(patches_1), 'r');hold on;end
% plot((1:size(stims,2))/fs,mean(stims)/piezo_scaling-piezo_offset, 'k');

% bottom = -25;
% for jj = 1:length(decs)-1
%     belongs = Spikes.spike_locs > incs(jj) & Spikes.spike_locs < decs(jj)+extra;
% if any(belongs)
%     plot([Spikes.spike_locs(belongs)-decs(jj)+2*extra Spikes.spike_locs(belongs)-decs(jj)+2*extra]/fs,...
%         [-jj+bottom -jj+0.5+bottom], 'k-');
% end
% end

% axis tight;% ylim([-60 -33]);
% xlim([0.495 0.52]);

if output_figure == 1
    try export_fig(fig2,[outputpath '/' filename(1:end-2) '.pdf'], '-pdf','-nocrop', '-r600' , '-painters', '-rgb');
    catch;end;
end
% datestr(now, 'mm_dd_yy')