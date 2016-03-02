function [spike_locs] = simple_spike_extractor_01_2015(data,unfiltered_data, piezo, sd_threshold)
%% simple version of a spike detector, requiring approval
spikeTemplateWidth = 50;
max_len = 400000;
if length(unfiltered_data) < max_len 
len = length(unfiltered_data)-100;
else
len = max_len -100;
end

%% find peaks for candidate locations  
if isempty(sd_threshold)
threshold_min = 0.5; threshold_max = 11;
thresholds = threshold_min:(threshold_max-threshold_min)/20:threshold_max;
else
thresholds = sd_threshold;
end

jj = 1;locs = [];
for ii = thresholds
locks = [];loccs = [];
% [pks,locks] = findpeaks(double(data),'minpeakheight',mean(data)+ii*std(data), 'minpeakdistance',floor(spikeTemplateWidth/2)-10);
[locks, pks] = peakfinder(double(data(1:len)),mean(data)+ii*std(data));%% slightly different algorithm;  [peakLoc] = peakfinder(x0,sel,thresh) returns the indicies of local maxima that are at least sel above surrounding data and larger (smaller) than thresh if you are finding maxima (minima).
 loccs = locks((locks> spikeTemplateWidth)); %%to prevent strange happenings, make sure that spikes do not occur right at the edges
locs(jj).spikes = loccs(loccs < (length(data)-spikeTemplateWidth));
jj = jj+1;
end

global Param;
if isempty(sd_threshold)
sliderGUI(data, unfiltered_data, piezo, locs, len, thresholds, []);
display('select a threshold and hit space');pause;
else 
Param = 1;
end



%% get all the spike locs using the correct std value
locks = [];loccs = [];spike_locs  = []; spike_history = [];
[locks, pks] = peakfinder(double(data),mean(data)+thresholds(Param)*std(data));%% slightly different algorithm;  [peakLoc] = peakfinder(x0,sel,thresh) returns the indicies of local maxima that are at least sel above surrounding data and larger (smaller) than thresh if you are finding maxima (minima).
 loccs = locks((locks> spikeTemplateWidth)); %%to prevent strange happenings, make sure that spikes do not occur right at the edges
spike_loc_candidates = loccs(loccs < (length(data)));
% [plocs, pinds] = find(piezo>0);
% only use points near the stimulus
% old_piezo = piezo;

incs = find(diff(piezo) > 3);
decs = find(diff(piezo) < -3);

pinds = [];extra = 1000;
for jj = 1:length(incs)
    pinds = [pinds incs(jj):decs(jj)+extra];
end

spike_loc_candidates = intersect(spike_loc_candidates,pinds);
plot_piezo = (piezo-mean(piezo))/max(piezo);

figure(65);clf; set(65, 'Position', [100 175 1000 1000],'color', 'w');
ind = 1;time_vec = (1:length(unfiltered_data))/10000;
for ii = 1:length(spike_loc_candidates)
    spike_flag = [];
    curr_candidate = spike_loc_candidates(ind);
    subplot(3,1,1);hold on;wind = 800;
    plot(time_vec(curr_candidate-wind:curr_candidate+wind),unfiltered_data(curr_candidate-wind:curr_candidate+wind));
    plot([time_vec(curr_candidate) time_vec(curr_candidate)] , [unfiltered_data(curr_candidate) unfiltered_data(curr_candidate)], 'ro');
    plot(time_vec(curr_candidate-wind:curr_candidate+wind), plot_piezo(curr_candidate-wind:curr_candidate+wind), 'r');
    plot(time_vec(spike_locs(spike_locs>curr_candidate-wind)), unfiltered_data(spike_locs(spike_locs>curr_candidate-wind)), 'go');axis tight;
    title({'press 0 if its a spike','press space if its not', 'press 1 to redo'});
    subplot(3,1,2);hold on;wind = 400;
    plot(time_vec(curr_candidate-wind:curr_candidate+wind),unfiltered_data(curr_candidate-wind:curr_candidate+wind));
    plot([time_vec(curr_candidate) time_vec(curr_candidate)] , [unfiltered_data(curr_candidate) unfiltered_data(curr_candidate)], 'ro');
    plot(time_vec(curr_candidate-wind:curr_candidate+wind), plot_piezo(curr_candidate-wind:curr_candidate+wind), 'r');
    plot(time_vec(spike_locs(spike_locs>curr_candidate-wind)), unfiltered_data(spike_locs(spike_locs>curr_candidate-wind)), 'go');axis tight;
    
    try 
    subplot(3,1,3);hold on;wind = 250;
    plot(time_vec(curr_candidate-wind:curr_candidate+wind),data(curr_candidate-wind:curr_candidate+wind));
    plot([time_vec(curr_candidate) time_vec(curr_candidate)] , [data(curr_candidate) data(curr_candidate)], 'ro');
    plot(time_vec(curr_candidate-wind:curr_candidate+wind), plot_piezo(curr_candidate-wind:curr_candidate+wind), 'r');
    plot(time_vec(spike_locs(spike_locs>curr_candidate-wind)), data(spike_locs(spike_locs>curr_candidate-wind)), 'go');axis tight;
    title(['candidate ' num2str(ii) 'of ' num2str(length(spike_loc_candidates))]);
    catch err;
    subplot(3,1,3);hold on;wind = 100;
    plot(time_vec(curr_candidate-wind:curr_candidate+wind),data(curr_candidate-wind:curr_candidate+wind));
    plot([time_vec(curr_candidate) time_vec(curr_candidate)] , [data(curr_candidate) unfiltered_data(curr_candidate)], 'ro');
    plot(time_vec(curr_candidate-wind:curr_candidate+wind), plot_piezo(curr_candidate-wind:curr_candidate+wind), 'r');
    plot(time_vec(spike_locs(spike_locs>curr_candidate-wind)), data(spike_locs(spike_locs>curr_candidate-wind)), 'go');axis tight;
    title(['candidate ' num2str(ii) 'of ' num2str(length(spike_loc_candidates))]);
    end

% spike_history = [0 1 1 1 1 1 1 1 1 1 1];spike_locs = [ 10:10:100];ind = 11;

get_key = waitforbuttonpress;
if get_key; spike_flag = get(gcf, 'CurrentCharacter');
    disp(spike_flag);end
if str2num(spike_flag) == 0
    spike_locs = [spike_locs; curr_candidate];
elseif str2num(spike_flag) == 1
    backprop = 5;
    if ind>backprop
     spike_locs(1+length(spike_locs)-sum(spike_history(end-(backprop-1):end)):end) = [];
%      sum(spike_history(end-(backprop-1):end))
     spike_history(end-(backprop-1):end) = [];
     ind = ind-backprop-1;
    else
    ind = 0;spike_locs = [];spike_history = [];
    end
   
end

if ~isscalar(str2num(spike_flag)) %% keep track of previous spike rankings
spike_history = [spike_history; 0];
elseif str2num(spike_flag) == 0 
spike_history = [spike_history; 1]; 
end

clf;
ind = ind+1;

end
end
