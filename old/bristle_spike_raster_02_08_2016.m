%% 02/08/2016 single trial analysis for bristle recording
clear all;
datapath = 'E:\ephys_data\';
spikepath = 'E:\ephys_data\spikes\';
outputpath = 'C:\Users\Tuthill Lab\Dropbox\tuthill lab shared\bristles\';

files = {'2016_02_16_0003.abf'};
ranges = {[]};

manual = 0;%% if spike sorting is to be done manually
output_figure = 0; %%make pdf figs?
other_piezo = 0; %% if 2nd piezo is the primary one, == 1; currently just one piezo
spike_detect = 1;%% 1 if new bristle spikes are to be detected, 0 if not
brist = 1;%% if brist = 1, plot bristle
patch = 0; %% if patch = 1, plot patch
redo = 0; %% to redo bristle spike detectiom
spikes.sd = 0.3; %%%made up standard deviation of the data (complicated by steady state deflection)

% for abf = 1: length(files) %% if averaging is to be done across files
abf = 1;
filename = files{abf};
[data,si,h]=abfload([datapath filename]);
fs = si/0.01;%% sampling rate
metadata = readtable(['E:\stim_data\' files{1}(1:end-4) '.txt'],'Delimiter','\t');
num_conditions = height(metadata);

spike_file = [spikepath filename(1:end-4) '_spikes.mat'];
patch_1 = data(:,1); %%patch primary signal
patch_2 = data(:,2); %%patch current signal
piezo_1 = -data(:,3);
piezo_1_sns = data(:,4);
bristle_1 = data(:,5);
cond_sig = data(:,6);

clear data;
if ~isempty(ranges{abf}); patch_1 = patch_1(ranges{abf});end
if ~isempty(ranges{abf}); patch_2 = patch_2(ranges{abf});end
if ~isempty(ranges{abf}); piezo_1 = piezo_1(ranges{abf});end
if ~isempty(ranges{abf}); piezo_1_sns = piezo_1_sns(ranges{abf});end
if ~isempty(ranges{abf}); bristle_1 = bristle_1(ranges{abf});end
patch_1_range = max(patch_1) - min(patch_1);

cutoff = 400;%%cutoff frequencies for high-pass filtering data
filts = cutoff/(fs/2);
[x,y] = butter(2,filts,'high');%%bandpass filter between 50 and 200 Hz
filtered_data = filter(x, y, bristle_1)';

if other_piezo == 1; piezo = piezo_2; else piezo = piezo_1; end
piezo(piezo < 0) = 0;%% assume all piezo movements are positive
if exist(spike_file) && redo ~=1;load(spike_file);spike_detect =0;end

time_vec = (1:length(patch_1))/fs;
secs = length(patch_1)/fs;
plot_range = 1*fs:secs*fs;

%% if necessary, pull out spikes
if spike_detect == 1
fly_start = 1;
baf = 0; %%number of samples to store before and after stimulus
spikes.present = 1; %%1 if its a spiking abf, empty if it's not
spikes.approval = 1; %%approval (== 0) to run without asking for the spike distance threshold
spikes.approval_time = 4;
spikes.threshold_factor = 1.4; %%threshold over which to find peaks (# of sds)
spikes.threshold = spikes.threshold_factor*spikes.sd;
spikes.spike_count_range = 1:length(bristle_1)+baf; %%range from which to extract spikes
spikes.filts = cutoff/(fs/2);
spikes.spikeDist_threshold = 10;
spikes.spikeTemplateWidth = 45; %%number of samples for the spike template
spikes.spikes2avg = 3;
spiketemplate = [];

Index = 1;
Spikes(Index).bristle_1    =   bristle_1';
Spikes(Index).piezo  =   piezo';

if manual == 1
[spike_locs] = simple_spike_extractor_02_08_2016(filtered_data(spikes.spike_count_range)', bristle_1(spikes.spike_count_range)',piezo(spikes.spike_count_range));
       Spikes(Index).num_spikes =   [Spikes(Index).num_spikes; length(spike_locs)];%%each row is a list of spike locations       
else
[spike_locs, spikesBelowThresh, spiketemplate, spikeDist_threshold] = extracellular_spike_extractor_02_08_2016(filtered_data(spikes.spike_count_range)', bristle_1(spikes.spike_count_range),...
   piezo(spikes.spike_count_range) ,spiketemplate, spikes.spikeTemplateWidth, spikes.spikes2avg, spikes.spikeDist_threshold, spikes.approval, spikes.approval_time);
end
Spikes(Index).num_spikes =   [length(spike_locs)];%%each row is a list of spike locations
       Spikes(Index).spike_ts =   [spikesBelowThresh'];%%each row is a single spike      
       Spikes(Index).spike_locs =   [spike_locs];%%the first column gives the spike location
       trial_num = size(Spikes(Index).bristle_1,1);
       Spikes(Index).spike_locs_inds =   [trial_num*ones(length(spike_locs),1)];%%the 2nd column gives the trial in which that spike location resides
       
       spikes.spikeDist_threshold = spikeDist_threshold; %%use the new threshold for the remainder of this abf
       bristle_spike_locs = zeros(1,length(bristle_1));bristle_spike_locs(Spikes(1).spike_locs) = 1; %% make a variable that shows all bristle spike locations

save([spikepath filename(1:end-4) '_spikes'], 'Spikes');
end

%% pull out each condition signal
expected_voltages = table2array(metadata(1:num_conditions,8));
cs_incs = find(diff(cond_sig) > 0.25);
for jj = 1:length(cs_incs)
    observed_trial_lengths(jj) = 100+find(diff(cond_sig(cs_incs(jj)+100:end)) < -0.25, 1);%%assuming that trial lasts more than 100 samples
    observed_voltages(jj) = round(10*mean(cond_sig(cs_incs(jj):cs_incs(jj)+observed_trial_lengths(jj))))/10;
    conditions(jj) = find(observed_voltages(jj) == expected_voltages);
end

for jj = 1:num_conditions  % create an array of the important variables for each condition
    trial_lengths(jj) = round(mean(observed_trial_lengths(conditions==jj)));%% trial lengths to use from here on out
    Data(abf,jj).patch_1 = [];
    Data(abf,jj).piezo_1 = [];    
    Data(abf,jj).bristle_1 = [];    
    Data(abf,jj).bristle_spikes_1 = [];
    Data(abf,jj).patch_spikes_1 = [];
    Data(abf,jj).start = [];
end

%% pull out each piezo stimulus
for cond = 1:length(conditions)
    
    start_points = cs_incs(conditions == cond);
    stop_points = start_points+trial_lengths(cond);
    incs = [];decs= [];this_piezo = [];this_bristle = [];this_patch = [];
    for rep = 1:sum(conditions == cond)
    
        this_piezo = piezo(start_points(rep):stop_points(rep));
        this_bristle = bristle_1(start_points(rep):stop_points(rep));
        this_patch = patch_1(start_points(rep):stop_points(rep));
        
        incs = find(diff(this_piezo) > 0.1);
        decs = find(diff(this_piezo) < -0.1);
        if decs(1) < incs(1);decs(1) = [];end
        if any(diff(incs) < 5); incs(diff(incs) < 10) = [];end %%make sure decs/incs are not too close
        if any(diff(decs) < 5); decs(diff(decs) < 10) = [];end
        if incs(length(incs)) > decs(length(decs)); incs(length(incs)) = [];end

%% if necessary, deal with light stimuli that have a faster duty cycle
        diff_incs = diff(incs);new_decs = [];new_incs = [];
        if any(diff_incs > 1 & diff_incs < 200) %% for cases where the light stimulus has a rapid duty cycle- use only the onset and offset
            inc_ups = find(diff_incs > 500);
            num_cycles = inc_ups(2)- inc_ups(1);
            new_incs = incs(1:num_cycles:end);
            for cyc = 1:length(new_incs)
                new_decs(cyc) = decs((num_cycles*cyc)-cyc);
                piezo(new_incs(cyc):new_decs(cyc)) = max(piezo);
            end
            decs = new_decs;
            incs = new_incs;
        end
   

%% parse into Data structure based on trials
        piezo_1_conds = [];
        extra = 1000;%% the extra data points to keep on either side of the trial
    for cycle = 1:(length(incs))    % search though all data and pull out relevant trials based on piezo onset
     current_start = incs(cycle);  % = the start position of each condition
     curr_Range = incs(cycle)-extra:incs(cycle)+(on_time)+2*extra; % creates a 1 x 2000 matrix (4 seconds at 500 hz)
     curr_piezo = abs(piezo_1(curr_Range));
     piezo_noise = std(curr_piezo(50:extra-50)); %%get the baseline noise to look for onset
     start = find(abs(curr_piezo) > piezo_noise*6, 1);%% set a threshold above which to look for onset, 6 sd's seems good
   
    cond_voltage(cycle) = mean(piezo(curr_Range(start+1:start+on_time-1)));
    srows = [];scols = []; svals = [];
        Index = cond(cycle);
        Data(abf,cond,rep).patch_1   =   [Data(abf,cond,rep).patch_1; patch_1(curr_Range(start - extra: start + on_time+extra))];
        Data(abf,cond,rep).piezo_1     =   [Data(abf,cond,rep).piezo_1; piezo_1(curr_Range(start- extra: start + on_time+extra))];
        Data(abf,cond,rep).bristle_1     =   [Data(abf,cond,rep).bristle_1; bristle_1(curr_Range(start- extra: start + on_time+extra))];
        Data(cond,rep).bristle_spikes_1   =   [Data(abf,cond,rep).bristle_spike_locs; bristle_spike_locs(curr_Range(start- extra: start + on_time+extra))];
        Data(cond,rep).start = start;
    end
    end
end
   
%% plot data
piezo_offset = -min(patch_1)+2;
piezo_scaling = 3;
add_to_bristle = 15;

fig2 = figure(2);clf;hold all;set(fig2,'color', 'w', 'position', [1 1 800 1050]);
cmap = colormap(autumn(length(decs)));
cmap2 = colormap(winter(length(decs)));

subplot(4,1,1); box off; hold all;title(filename, 'interpreter', 'none');
p(1,:) = plot(time_vec, piezo_1/piezo_scaling-piezo_offset, 'b', 'linewidth', 1);
if patch ==1;plot(time_vec, patch_1, 'k', 'linewidth', 1);end
if brist ==1; plot(time_vec, bristle_1+add_to_bristle, 'linewidth', 1, 'color', [0.3 0.3 0.3]);end
% plot(time_vec, cond_sig)

subplot(4,1,2); hold all;box off;title('all trials');
buffer = 200;
for jj = 1:length(decs)
ext = time_vec(1:decs(jj)-incs(jj)+2*buffer+1);
 if brist == 1; plot(ext,add_to_bristle+bristle_1((incs(jj)-buffer):(decs(jj)+buffer)), 'linewidth', 0.75, 'color', cmap2(jj,:));end
if patch == 1; plot(ext, patch_1((incs(jj)-200):(decs(jj)+200)), 'linewidth', 0.75, 'color', cmap(jj,:));end
plot(ext, piezo((incs(jj)-buffer):(decs(jj)+buffer))/piezo_scaling-piezo_offset, 'k', 'linewidth', 0.75);
end

bottom = -25+add_to_bristle;
for jj = 1:length(decs)-1
    belongs = Spikes.spike_locs > incs(jj) & Spikes.spike_locs < decs(jj)+extra;
if any(belongs)
    plot([Spikes.spike_locs(belongs)-decs(jj)+extra+buffer Spikes.spike_locs(belongs)-decs(jj)+extra+buffer]/fs,...
        [-jj+bottom -jj+0.5+bottom], 'k-');
end
end

axis tight;% ylim([-60 -33]);
% xlim([0 0.55]);


% subplot(4,1,3:4);hold all;box off;title('averages');
% if brist ==1; confplot((1:size(bristles_1,2))/fs,mean(bristles_1)-10, std(bristles_1), std(bristles_1), 'b');hold on;end
% if patch ==1; confplot((1:size(patches_1,2))/fs,mean(patches_1), std(patches_1), std(patches_1), 'r');hold on;end
% plot((1:size(stims,2))/fs,mean(stims)/piezo_scaling-piezo_offset, 'k');

subplot(4,1,3:4);hold all;box off;title('averages');
if brist ==1; plot((1:size(bristles_1,2))/fs,mean(bristles_1)-10, 'b');hold on;end
if patch ==1; plot((1:size(patches_1,2))/fs,mean(patches_1), 'r');hold on;end
plot((1:size(stims,2))/fs,mean(stims)/piezo_scaling-piezo_offset, 'k');

bottom = -35;
for jj = 1:length(decs)
    belongs = Spikes.spike_locs > incs(jj) & Spikes.spike_locs < decs(jj)+extra;
        total_spikes(jj) = sum(belongs);

if any(belongs)
    plot([Spikes.spike_locs(belongs)-decs(jj)+2*extra Spikes.spike_locs(belongs)-decs(jj)+2*extra]/fs,...
        [-jj+bottom -jj+0.5+bottom], 'k-');
end
end

% axis tight;% ylim([-60 -33]);
% xlim([0.495 0.52]);

legend('light', 'mechanical', 'location', 'best');


fig9 = figure(9);clf;hold all;set(fig9,'color', 'w', 'position', [1 1 800 1050]);
buffer = 1000;
for jj = 1:length(decs)
ext = time_vec(1:decs(jj)-incs(jj)+2*buffer+1);
 if brist == 1; plot(ext,add_to_bristle+bristle_1((incs(jj)-buffer):(decs(jj)+buffer))-jj*5, 'linewidth', 0.75, 'color', cmap2(jj,:));end
plot(ext, piezo((incs(jj)-buffer):(decs(jj)+buffer))/piezo_scaling-piezo_offset, 'k', 'linewidth', 0.75);
end

    export_fig(fig9,[outputpath '/' filename(1:end-2) '_exam.pdf'], '-pdf','-nocrop', '-r600' , '-painters', '-rgb');



if output_figure == 1
    try export_fig(fig2,[outputpath '/' filename(1:end-2) '.pdf'], '-pdf','-nocrop', '-r600' , '-painters', '-rgb');
    export_fig(fig4,[outputpath '/' filename(1:end-2) '_raster.pdf'], '-pdf','-nocrop', '-r600' , '-painters', '-rgb');
    export_fig(fig5,[outputpath '/' filename(1:end-2) '_means.pdf'], '-pdf','-nocrop', '-r600' , '-painters', '-rgb');

    catch;end;
end
% datestr(now, 'mm_dd_yy')