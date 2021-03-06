%% 02/08/2016 single trial analysis for bristle recording
clear all;

%% init file and paths
files = {'2016_03_02_0006.abf'};
ranges = {[1:3.8e6],[],[],[],[],[],[],[]};

rootfolder = 'C:\Data\';
datapath = [rootfolder 'ephys_data\'];
spikepath = [rootfolder 'spikes\'];
stimpath = [rootfolder 'stim_metadata\'];
summpath = [rootfolder 'ephys_metadata\'];
outputpath = [rootfolder 'figures\'];

%% init analysis params
manual = 0;%% if spike sorting is to be done manually
output_figure = 0; %%make pdf figs?
other_piezo = 0; %% if 2nd piezo is the primary one, == 1; currently just one piezo
spike_detect = 1;%% 1 if new bristle spikes are to be detected, 0 if not
brist = 1;%% if brist = 1, plot bristle
patch = 0; %% if patch = 1, plot patch
redo_spikes = 0; %% to redo_spikes bristle spike detectiom
spikes.sd = 0.3; %%%made up standard deviation of the data (complicated by steady state deflection)
firing_bin_ms = 100;%% window in which to calculate firing rate (in samples)

for abf = 1: length(files) %% if averaging is to be done across files
filename = files{abf};
metadata = readtable([stimpath files{1}(1:end-4) '.txt'],'Delimiter','\t');
num_conditions = height(metadata);
num_cycles = mode(table2array(metadata(1:num_conditions,7)));
num_reps = mode(table2array(metadata(1:num_conditions,5)));
expected_voltages = table2array(metadata(1:num_conditions,8));

%% import data
[data,si,h]=abfload([datapath filename]);
spike_file = [spikepath filename(1:end-4) '_spikes.mat'];
patch_1 = data(:,1); %%patch primary signal
patch_2 = data(:,2); %%patch current signal
piezo_1 = -data(:,3);
piezo_1_sns = data(:,4);
bristle_1 = data(:,5);
cond_sig = data(:,6);
fs = si/0.01;%% sampling rate
firing_bin = firing_bin_ms*(fs/1000); 

clear data;
if ~isempty(ranges{abf}); patch_1 = patch_1(ranges{abf});end
if ~isempty(ranges{abf}); patch_2 = patch_2(ranges{abf});end
if ~isempty(ranges{abf}); piezo_1 = piezo_1(ranges{abf});end
if ~isempty(ranges{abf}); piezo_1_sns = piezo_1_sns(ranges{abf});end
if ~isempty(ranges{abf}); bristle_1 = bristle_1(ranges{abf});end
if ~isempty(ranges{abf}); cond_sig = cond_sig(ranges{abf});end
patch_1_range = max(patch_1) - min(patch_1);

%% filter bristle data
% cutoff_low = 100;
% cutoff_high = 1000;%%cutoff frequencies for filtering bristle recording data
% filts = [cutoff_low cutoff_high]/(fs/2);
% [x,y] = butter(2,filts,'stop');%%bandpass filter between 50 and 200 Hz

cutoff_high = 100;%%cutoff frequencies for filtering bristle recording data
filts = cutoff_high/(fs/2);
[x,y] = butter(4,filts,'high');%%bandpass filter between 50 and 200 Hz
filtered_data_high = filter(x, y, bristle_1)';

cutoff_low = 600;%%cutoff frequencies for filtering bristle recording data
filts = cutoff_low/(fs/2);
[x,y] = butter(2,filts,'low');%%bandpass filter between 50 and 200 Hz
filtered_data = filter(x, y, filtered_data_high);

if other_piezo == 1; piezo = piezo_2; else piezo = piezo_1; end
piezo(piezo < 0) = 0;%% assume all piezo movements are positive
if exist(spike_file) && redo_spikes ~=1;load(spike_file);spike_detect =0;end

%% if necessary, pull out spikes
if spike_detect == 1
fly_start = 1;
baf = 0; %%number of samples to store before and after stimulus
spikes.present = 1; %%1 if its a spiking abf, empty if it's not
spikes.approval = 1; %%approval (== 0) to run without asking for the spike distance threshold
spikes.approval_time = 10;
spikes.threshold_factor = 1.4; %%threshold over which to find peaks (# of sds)
spikes.threshold = spikes.threshold_factor*spikes.sd;
spikes.spike_count_range = 1:length(bristle_1)+baf; %%range from which to extract spikes
spikes.spikeDist_threshold = 10;
spikes.spikeTemplateWidth = 50; %%number of samples for the spike template
spikes.spikes2avg = 3;
spiketemplate = [];

Index = 1;
Spikes(Index).bristle_1    =   bristle_1';
Spikes(Index).piezo  =   piezo';

if manual == 1
[spike_locs] = simple_spike_extractor_02_08_2016(filtered_data(spikes.spike_count_range)', bristle_1(spikes.spike_count_range)',piezo(spikes.spike_count_range));
       Spikes(Index).num_spikes =   [Spikes(Index).num_spikes; length(spike_locs)];%%each row is a list of spike locations       
else
[spike_locs, spikesBelowThresh, spiketemplate, spikeDist_threshold, vars] = extracellular_spike_extractor_03_02_2016(bristle_1(spikes.spike_count_range),...
   piezo(spikes.spike_count_range) ,spiketemplate, spikes.spikeTemplateWidth, spikes.spikes2avg, spikes.spikeDist_threshold, spikes.approval, spikes.approval_time);
end
Spikes(Index).num_spikes =   [length(spike_locs)];%%each row is a list of spike locations
       Spikes(Index).spike_ts =   [spikesBelowThresh'];%%each row is a single spike      
       Spikes(Index).spike_locs =   [spike_locs];%%the first column gives the spike location
       trial_num = size(Spikes(Index).bristle_1,1);
       Spikes(Index).spike_locs_inds =   [trial_num*ones(length(spike_locs),1)];%%the 2nd column gives the trial in which that spike location resides
       
       Spikes.spikeDist_threshold = spikeDist_threshold; %%use the new threshold for the remainder of this abf
       Spikes.bristle_spike_locs = zeros(1,length(bristle_1));
       Spikes.bristle_spike_locs(Spikes(1).spike_locs) = 1; %% make a variable that shows all bristle spike locations
       Spikes.vars = vars;%% filtering parameters
save([spikepath filename(1:end-4) '_spikes'], 'Spikes');
end

%% pull out each condition signal
cond_sig(cond_sig<0) = 0;%% in case transients dip below zero
cs_incs = find(diff(cond_sig) > 0.25);

for jj = 1:length(cs_incs)
    observed_trial_lengths(jj) = 100+find(diff(cond_sig(cs_incs(jj)+100:end)) < -0.25, 1);%%assuming that trial lasts more than 100 samples
    observed_voltages(jj) = round(10*mean(cond_sig(cs_incs(jj):cs_incs(jj)+observed_trial_lengths(jj))))/10;
    conditions(jj) = find(observed_voltages(jj) == expected_voltages);
end

for cond = 1:num_conditions  % create an array of the important variables for each condition
       trial_lengths(cond) = round(mean(observed_trial_lengths(conditions==cond)));%% trial lengths to use from here on out
    for cycle = 1:num_cycles
    Data(abf,cond,cycle).patch_1 = [];
    Data(abf,cond,cycle).piezo_1 = [];    
    Data(abf,cond,cycle).bristle_1 = [];    
    Data(abf,cond,cycle).bristle_spikes_1 = [];
    Data(abf,cond,cycle).bristle_spike_rate = [];
    Data(abf,cond,cycle).patch_spikes_1 = [];
    Data(abf,cond,cycle).start = [];
    end
end

%% pull out and parcel each piezo stimulus
for cond = 1:num_conditions
    start_points = cs_incs(conditions == cond);
    stop_points = start_points+trial_lengths(cond);
    incs = [];decs= [];thisrep_piezo = [];thisrep_bristle = [];thisrep_patch = [];
    
    for rep = 1:length(start_points)
        
        thisrep_start = start_points(rep);thisrep_stop = stop_points(rep);
        thisrep_piezo = piezo(thisrep_start:thisrep_stop);
        thisrep_bristle = bristle_1(thisrep_start:thisrep_stop);
        thisrep_patch = patch_1(thisrep_start:thisrep_stop);
        
        incs = find(diff(thisrep_piezo) > 0.1);
        thisrep_piezo(incs+1) = thisrep_piezo(incs);%% a hack to deal with transients
        decs = find(diff(thisrep_piezo) < -0.1);
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
   
    on_time = mode(decs-incs);
    off_time = mode(incs(2:end)-decs(1:end-1));
    
%% parse into Data structure based on trials
        extra = 1000;%% the extra data points to keep on either side of the trial
    for cycle = 1:(length(incs))    % search though all data and pull out relevant trials based on piezo onset
     these_bristle_spikes = [];
     curr_Range = incs(cycle):incs(cycle)+(on_time); % creates a 1 x 2000 matrix (4 seconds at 500 hz)
     cycle_start = thisrep_start+incs(cycle);
%      curr_piezo = abs(piezo_1(cycle_start:cycle_start+on_time));
%      piezo_noise = std(curr_piezo(50:extra-50)); %%get the baseline noise to look for onset

        Data(abf,cond,cycle).patch_1   =   [Data(abf,cond,cycle).patch_1; patch_1(cycle_start - extra: cycle_start + on_time+extra)'];
        Data(abf,cond,cycle).piezo_1     =   [Data(abf,cond,cycle).piezo_1; piezo_1(cycle_start- extra: cycle_start + on_time+extra)'];
        Data(abf,cond,cycle).bristle_1     =   [Data(abf,cond,cycle).bristle_1; bristle_1(cycle_start- extra: cycle_start + on_time+extra)'];
        
        these_bristle_spikes = Spikes.bristle_spike_locs(cycle_start- extra:cycle_start + on_time+extra);
        Data(abf, cond, cycle).bristle_spikes_1   =   [Data(abf,cond,cycle).bristle_spikes_1; these_bristle_spikes'];
        Data(abf,cond, cycle).start = cycle_start;
        
        smoothed_rate = zeros(1,length(Data(abf,cond,cycle).patch_1));
        for discs = 1+firing_bin/2:length(Data(abf,cond,cycle).patch_1)-firing_bin/2
        smoothed_rate(1,discs) = sum(these_bristle_spikes(discs-firing_bin/2:discs+firing_bin/2))*(fs/firing_bin);
        end
        Data(abf, cond, cycle).bristle_spike_rate = [Data(abf, cond, cycle).bristle_spike_rate; smoothed_rate];
        
    end
    end
end
   
save([summpath filename(1:end-4) '_summ'], 'Data');
end