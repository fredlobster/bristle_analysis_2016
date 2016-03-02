function [spike_locs] = simple_spike_extractor_07_2014(data,unfiltered_data, piezo)
%% simple version of a spike detector, requiring approval
spikeTemplateWidth = 60;
max_len = 400000;
if length(unfiltered_data) < max_len 
len = length(unfiltered_data)-1000;
else
len = max_len -1000;
end

%% find peaks for candidate locations  
%% 9_23_2012, invert the signal to look for troughs- better for extracellular spikes?
threshold_min = 3; threshold_max = 11;
thresholds = threshold_min:(threshold_max-threshold_min)/20:threshold_max;
jj = 1;locs = [];
for ii = thresholds
locks = [];loccs = [];
% [pks,locks] = findpeaks(double(data),'minpeakheight',mean(data)+ii*std(data), 'minpeakdistance',floor(spikeTemplateWidth/2)-10);
[locks, pks] = peakfinder(double(data(1:len)),mean(data)+ii*std(data));%% slightly different algorithm;  [peakLoc] = peakfinder(x0,sel,thresh) returns the indicies of local maxima that are at least sel above surrounding data and larger (smaller) than thresh if you are finding maxima (minima).
 loccs = locks((locks> spikeTemplateWidth)); %%to prevent strange happenings, make sure that spikes do not occur right at the edges
locs(jj).spikes = loccs(loccs < (length(data)-spikeTemplateWidth));
jj = jj+1;
end

if isempty(locs)%%make sure there are points that exceed the threshold
display('no points above threshold'); detected_spike_locs = []; spikesBelowThresh = []; spikeTemplate = template;
figure(9); set(9, 'Position', [200 175 800 500],'color', 'w');
subplot(3,1,1);plot(unfiltered_data); title('unfiltered data'); hold on;
subplot(3,1,2:3);plot(data); title('No spikes??'); pause(1);
return   
end

global Param;
% sliderGUI(data, unfiltered_data, piezo, locs, len, thresholds);
sliderGUI(data, unfiltered_data, piezo, locs, len, thresholds, []);

display('select a threshold and hit space');pause;

%% get all the spike locs using the correct std value
locks = [];loccs = [];spike_locs  = [];
[locks, pks] = peakfinder(double(data),mean(data)+thresholds(Param)*std(data));%% slightly different algorithm;  [peakLoc] = peakfinder(x0,sel,thresh) returns the indicies of local maxima that are at least sel above surrounding data and larger (smaller) than thresh if you are finding maxima (minima).
 loccs = locks((locks> spikeTemplateWidth)); %%to prevent strange happenings, make sure that spikes do not occur right at the edges
spike_loc_candidates = loccs(loccs < (length(data)-1000));
plot_piezo = (piezo-mean(piezo))/max(piezo);
[plocs, pinds] = find(piezo>0);
spike_loc_candidates = intersect(spike_loc_candidates,pinds);

figure(65); set(65, 'Position', [200 175 600 1000],'color', 'w');
for ii = 1:length(spike_loc_candidates)
    spike_flag = [];
    curr_candidate = spike_loc_candidates(ii);
    subplot(4,1,1);hold on;wind = 1000;
    plot(unfiltered_data(curr_candidate-wind:curr_candidate+wind));
    plot([wind wind] , [unfiltered_data(curr_candidate) unfiltered_data(curr_candidate)], 'ro');
    plot(plot_piezo(curr_candidate-wind:curr_candidate+wind), 'r');
    title('press 0 if its a spike');
    subplot(4,1,2);hold on;wind = 600;
    plot(unfiltered_data(curr_candidate-wind:curr_candidate+wind));
    plot([wind wind] , [unfiltered_data(curr_candidate) unfiltered_data(curr_candidate)], 'ro');
    plot(plot_piezo(curr_candidate-wind:curr_candidate+wind), 'r');
    
    subplot(4,1,3);hold on;title('filtered')
    try 
    lead = 600;
    plot(data(curr_candidate-lead:curr_candidate+lead));
    plot([lead lead] , [data(curr_candidate) data(curr_candidate)], 'ro');
    plot(plot_piezo(curr_candidate-lead:curr_candidate+lead), 'r');
    catch err;
    lead = 300;
    plot(data(curr_candidate-lead:curr_candidate+lead));
    plot([lead lead] , [data(curr_candidate) data(curr_candidate)], 'ro');
    plot(plot_piezo(curr_candidate-lead:curr_candidate+lead), 'r');
    end
%    subplot(4,1,4);hold on;
%     plot(unfiltered_data(curr_candidate-500:curr_candidate+500));
%     plot([1000 1000] , [unfiltered_data(curr_candidate) unfiltered_data(curr_candidate)], 'ro');

get_key = waitforbuttonpress;
if get_key; spike_flag = get(gcf, 'CurrentCharacter');
    display(spike_flag);end
if str2num(spike_flag) == 0
    spike_locs = [spike_locs; curr_candidate];
end
clf;
end

end
