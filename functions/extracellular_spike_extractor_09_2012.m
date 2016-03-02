function [detected_spike_locs, spikesBelowThresh, spikeTemplate, spikeDistThreshold] = extracellular_spike_extractor_09_2012(data,...
    unfiltered_data, threshold, template, spikeTemplateWidth,  numtemps, Distance_threshold, approval, approval_time)
%% Code to match a given spike waveform across a recording (1-D); basic idea is to estimate warping distance via DTW and convert the  estimated warp distance into a probability via a kernel
% John Tuthill, Parvez Ahammad, JFRC/HHMI; Last updated: 04/28/2011
detected_spike_locs = []; spikesBelowThresh = [];
spikeDistThreshold  = Distance_threshold;

%% find peaks for candidate locations  
%% 9_23_2012, invert the signal to look for troughs- better for extracellular spikes?
[pks,locks] = findpeaks(-double(data),'minpeakheight',mean(data)+threshold, 'minpeakdistance',floor(spikeTemplateWidth/3));

loccs = locks((locks > spikeTemplateWidth)); %%to prevent strange happenings, make sure that spikes do not occur right at the edges
locs = loccs(loccs < (length(data)-spikeTemplateWidth))

if isempty(locs)%%make sure there are points that exceed the threshold
display('no points above threshold'); detected_spike_locs = []; spikesBelowThresh = []; spikeTemplate = template;
figure(9); set(9, 'Position', [200 175 800 500],'color', 'w');
subplot(3,1,1);plot(unfiltered_data); title('unfiltered data'); hold on;
subplot(3,1,2:3);plot(data); title('No spikes??'); pause(1);
return   
end

if ~isempty(template)
spikeTemplate = template;
spikeTemplateWidth = length(spikeTemplate);
else
    
template_center = [];
fig = figure(9); set(9, 'Position', [200 100 800 600],'color', 'w');
subplot(3,1,1);plot(unfiltered_data); title('unfiltered data'); hold on;
subplot(3,1,2:3);plot(data); title('Click on a spike and then press enter'); hold on;
plot((locs), zeros(1,length(locs)),'ro');
xlabel(['sd = ' num2str(std(data))]);

dcm_obj = datacursormode(fig);
set(dcm_obj,'DisplayStyle','datatip','SnapToDataVertex','off','Enable','on')

spikes_flag = [];
while isempty(spikes_flag)          % Wait while the user does this.
   for hh = 1:numtemps
   spikes_flag = input('\n If there are no spikes, press 1...');
    try
    figure(9);template_center = [];
    c_info = getCursorInfo(dcm_obj);
    template_center = c_info.DataIndex;
    spikeTemplate(hh,:) = data(template_center-(0.5*spikeTemplateWidth):template_center+(0.5*spikeTemplateWidth));
    catch err
    display('no spikes selected'); detected_spike_locs = []; spikesBelowThresh = []; spikeTemplate = template;
   end
end
    spikes_flag = 0;
        if size(spikeTemplate,1)>1; spikeTemplate = mean(spikeTemplate,1);end
end

if isempty(spikeTemplate); detected_spike_locs = []; spikesBelowThresh = []; 
    return; end;

close;
figure(12), plot(spikeTemplate), title('Template Waveform: wait one second');pause(1);close;
end

spikeTemplateWidth = length(spikeTemplate);

%% pool the detected spike candidates and do template matching
targetSpikeDist = zeros(size(locs));
counter = 1;
for i=1:length(locs)

    if min(locs(i)+spikeTemplateWidth/2,length(data)) - max(locs(i)-spikeTemplateWidth/2,0)< spikeTemplateWidth
        continue
    else
%         curSpikeTarget = data(max(locs(i)-floor(spikeTemplateWidth/2),0)+1:min(locs(i)+floor(spikeTemplateWidth/2),length(data)),1);  %% old method
%         detectedSpikeCandidates(:,counter) = data(max(locs(i)-floor(spikeTemplateWidth/2),0): min(locs(i)+floor(spikeTemplateWidth/2),length(data)),1);
        curSpikeTarget = data(max(locs(i)-floor(spikeTemplateWidth/2),0)+1: min(locs(i)+floor(spikeTemplateWidth/2),length(data)));
        detectedSpikeCandidates(:,counter) = data(max(locs(i)-floor(spikeTemplateWidth/2),0): min(locs(i)+floor(spikeTemplateWidth/2),length(data)));
        norm_curSpikeTarget = curSpikeTarget/max(curSpikeTarget);norm_spikeTemplate = spikeTemplate/max(spikeTemplate);
        % figure(34); hold all; plot(curSpikeTarget, 'b');plot(norm_curSpikeTarget, 'r');pause;
        [targetSpikeDist(i), foo,bar] = dtw_WarpingDistance(norm_curSpikeTarget, norm_spikeTemplate);
        clear foo bar;      
    end
    counter = counter+1;
end

% estimate spike probabilities at candidate locations
spikeProbs = zeros(size(locs));
for i=1:length(locs)
   if min(locs(i)+spikeTemplateWidth/2,length(data)) - max(locs(i)-spikeTemplateWidth/2,0)< spikeTemplateWidth
       continue
   else
       spikeProbs(i) = exp( -(abs(targetSpikeDist(i)-mean(targetSpikeDist))) / (2*var(targetSpikeDist)) );
   end
end

%%
figure(11); clf; set(11, 'Position', [0 0 1000 700],'color', 'w');
accept_spikeDistThreshold_flag = 0;
while accept_spikeDistThreshold_flag == 0
detected_spike_locs = [];clear spikesAboveThresh spikesBelowThresh spikesAboveThresh_aligned spikesBelowThresh_aligned;
figure(11); clf;
subplot(3,2,1); plot(detectedSpikeCandidates), title('All detected spike candidates');hold on; plot(spikeTemplate, 'k', 'linewidth', 4)
subplot(3,2,3), hist(targetSpikeDist,20), title('Histogram of distances')
subplot(3,2,5);plot(unfiltered_data-mean(unfiltered_data)), hold on; plot(detected_spike_locs, zeros(1,length(detected_spike_locs)),'ro');axis off;

    spikesAboveThresh = []; spikesBelowThresh = [];spikeDistThreshold = [];detected_spike_locs = [];
    
    if isempty(Distance_threshold)
    subplot(3,2,5);title('Enter the new spike distance threshold:');hold all;
    spikeDistThreshold = input('\n \n Enter the spike distance threshold:');
    else
    spikeDistThreshold  = Distance_threshold;
    end
    i =1;
    aboveThreshCounter =1;
    belowThreshCounter=1;
    
    while i<= length(locs)
                if length(data(max(locs(i)-floor(spikeTemplateWidth/2))+1: min(locs(i)+floor(spikeTemplateWidth/2),length(data)))) >= spikeTemplateWidth-1
            curSpikeTarget = data(max(locs(i)-floor(spikeTemplateWidth/2))+1: min(locs(i)+floor(spikeTemplateWidth/2),length(data)));
                    %         if length(data(max(locs(i)-floor(spikeTemplateWidth/2),0)+1: min(locs(i)+floor(spikeTemplateWidth/2),length(data)),1)) >= spikeTemplateWidth-1
%             curSpikeTarget = data(max(locs(i)-floor(spikeTemplateWidth/2),0)+1: min(locs(i)+floor(spikeTemplateWidth/2),length(data)),1);
            if targetSpikeDist(i)> spikeDistThreshold
                spikesAboveThresh(:, aboveThreshCounter)= curSpikeTarget;
                aboveThreshCounter = aboveThreshCounter+1;
            else
                spikesBelowThresh(:, belowThreshCounter)= curSpikeTarget;
                belowThreshCounter = belowThreshCounter+1;
                detected_spike_locs = [detected_spike_locs;locs(i)];
            end
        end
        i=i+1;
    end
    % update the counters
    belowThreshCounter = belowThreshCounter-1;
    aboveThreshCounter = aboveThreshCounter-1;
    
%     ind = detected_spikes<spikeDistThreshold;
%     finalSpikeLocs = detected_spike_locs;
    
    figure(11);
    subplot(3,2,2), plot(spikesAboveThresh), title(sprintf('Spikes above %f: number = %d', spikeDistThreshold, size(spikesAboveThresh,2)))
    subplot(3,2,4), plot(spikesBelowThresh), title(sprintf('Spikes below %f: number = %d', spikeDistThreshold, size(spikesBelowThresh,2)))
    subplot(3,2,6), plot(data-mean(data)), hold on; plot(detected_spike_locs, zeros(1,length(detected_spike_locs)),'ro'); title ('Estimated spikes');
    if approval== 0
      accept_spikeDistThreshold_flag = 1;
    elseif approval == 1
       subplot(3,2,5), axis off; title({'Are you happy with the spike output?' 'Hit 0 if No'}); hold all;
       plot(unfiltered_data-mean(unfiltered_data)), hold on; plot(detected_spike_locs, zeros(1,length(detected_spike_locs)),'ro');
       if ~exist('approval_time','var')
           approval_time = 1;
       end
       accept_spikeDistThreshold_flag = timeinput(approval_time); %%returns 1 if no input within 1 second
       Distance_threshold = [];
    elseif approval == 2
    subplot(3,2,5), axis off; title({'Are you happy with the spike output?' 'Hit Enter for yes or press 1 for No'}); hold all;
    plot(unfiltered_data-mean(unfiltered_data)), hold on; plot(detected_spike_locs, zeros(1,length(detected_spike_locs)),'ro');
    accept_spikeDistThreshold_flag = input('\n Are you happy with the spike output? \n Hit Enter for yes or press 0 for No:');
    Distance_threshold = [];
    end
end
end
