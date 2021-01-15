%% Analysis script to analyze single subjects of the BeeG dataset

function do_singlesubject_analysis(sub, subjectnumber, subjectlist, data_dir)

%% Input of the dataset

% Start by running the subject specific details script to find relevant
% info

details_script_name = erase(sub, 'sub-');
details_script_name = ['details_sub' details_script_name];
run(details_script_name)

if ~exist(output_dir, 'dir')
   mkdir(output_dir);    
end

%% Find the data and headerfiles

cfg.dataset = [data_dir filesep sub filesep 'eeg' filesep sub '_task-audiovisual_eeg.vhdr' ];
hdr         = ft_read_header(cfg.dataset);

%% Now we can extract the trials from the events.tsv

t              = readtable([data_dir filesep sub filesep 'eeg' filesep sub '_task-audiovisual_events.tsv'], 'FileType', 'text');
trl            = t( : , 3:8);

% The variable raw now contains the eeg tsv info:
% Column 1: begsample
% Column 2: endsample
% Column 3: offset
% Column 4: marker number
% Column 5: stimulus
% Column 6: location of the bee

% Then we re-define trials. We are only interested in bees and cues
% And take a time window of -500 ms and + 1000 ms

stimuli             = {'bee', 'update-cue', 'no-update-cue'};
str                 = string(trl{:,5});
index_stimuli       = ismember(str, stimuli);
trl_cueonly         = trl(index_stimuli, :); % now only bee and cue trials are left

pre_stim_samples    = round(0.5 * hdr.Fs); % The pre-stim period is 0.5 s
post_stim_samples   = round(1 * hdr.Fs);   % The post-stim period is 1 s

% If the first trial starts to early, we remove it
if trl_cueonly{1,1} < pre_stim_samples
    trl_cueonly     = trl_cueonly(2:end, :);
    trl_new         = trl_cueonly;
else
    trl_new         = trl_cueonly;
end

trl_new{:,1}        = trl_cueonly{:,1}-pre_stim_samples; % We extract the pre stim samples from begsample (the first colum of trl_cueonly) to find the new begsample
trl_new{:,2}        = trl_cueonly{:,1}+post_stim_samples; % We add the post stim samples to begsample to find the new endsample
trl_new{:,3}        = - pre_stim_samples;

% Then add the new trials to cfg
cfg.trl             = trl_new;

%% Pre-processing on the episodes that are defined as trials

% We perform pre-processing similarly to the published work on this study:
% Kayhan et al. Developmental Cognitive Neuroscience, 2019

% Concretely: 5 sec padding for high pass filtering at 1 Hz and baseline
% correction on the entire window

% Set channels
cfg.channel             = ft_channelselection('EEG', hdr.label); % Only the EEG channels

% Baseline-correction options
cfg.bpfilter            = 'yes';  % band-pass filtering  
cfg.bpfreq              = [1 30]; % filter between 1-30 Hz
cfg.padding             = 5;
cfg.demean              = 'yes';
data                    = ft_preprocessing(cfg);

%% Artefact rejection

% We start the artefact rejection, consisting of three phases:
% 1) A rough visual rejection
% 2) ICA
% 3) Finetuning by rejecting trials via visual artefact rejection

% This is only necessary once, if it has been done, we skip this step

if do_artefact_rejection==1
    % we do the artefact rejection

    %% Artefact rejection part 1: a rough visual data rejection using ft_rejectvisual

    % We start with trial view to get an overview of the data, then summary
    % view to reject channels or trials that are way off

    cfg                 = [];
    cfg.method          = 'trial';  % Or switch to summary if needed
    cfg.keepchannel     = 'nan';    % when rejecting channels, values are replaced by NaN
    cfg.alim            = [-100, 100];   
    data_artefact_1     =  ft_rejectvisual(cfg, data);
    
    % Let's find the bad channels so we can interpolate them later
    badchannels         = find(all(isnan(data_artefact_1.trial{1}), 2));

    % Then we find the rejected trials
    badtrial_times      = data_artefact_1.cfg.artfctdef.trial.artifact; % this matrix contains start and end times of each rejected trial

    %% Artefact rejection part 2: ICA

    % We do this only on the 'good channels, so first we remove the bad channels'

    channels                 = data_artefact_1.label;
    channels(badchannels, :) = [];

        % Then we set up the cfg struct and perform ICA

    cfg                      = [];
    cfg.channel              = channels;
    cfg.method               = 'runica'; % this is the default
    comp                     = ft_componentanalysis(cfg, data_artefact_1);

    % Then we plot the first 20 components
    cfg                      = [];
    cfg.component            = 1:20; % We plot the first 20 components      
    cfg.layout               = 'EEG1010.lay';
    cfg.comment              = 'no';
    cfg.zlim                 = 'maxabs'; 
 
    ft_topoplotIC(cfg, comp)
    
    disp('')
    input('Press Enter to continue')
    close all
     
    %In the end, plot the time course of all components again to check whether the components 
    % you want to remove indeed look like artifacts over time
    
    cfg                      = [];
    cfg.layout               = 'EEG1010.lay'; % specify the layout file that should be used for plotting
    cfg.viewmode             = 'component';
    ft_databrowser(cfg, comp)
    
    disp('')
    input('Press Enter to continue')
    close all
 
    % Remove the bad components and backproject the data  
    prompt                       = 'Which components do you want to reject [enter row vector]';
    rejcom                       = input(prompt);

    % Now we exlude all channels excluded during the first visual artefact rejection (containing NaNs)

    exclude                      = find(~ismember(data_artefact_1.label, channels));
    tempdata                     = data_artefact_1;
    for trl=1:length(tempdata.trial)
        tempdata.trial{1,trl}(exclude,:) = [];
        tempdata.label = channels;
    end

    cfg                          = [];
    cfg.component                = rejcom; % to be removed components
    cfg.channel                  = channels;
    tempcleandata                = ft_rejectcomponent(cfg, comp, tempdata);

    % Then we copy these back into the full data structure with the 'bad
    % channels' as nans

    cleandata                    = data_artefact_1;
    for trl = 1:length(cleandata.trial) % Loop through trials
        in=1; 
        for ch = 1:size(cleandata.trial{1,trl},1) % Loop throuch channels
            if ~isnan(cleandata.trial{1,trl}(ch,1)) % if the channels is not a bad channels with only isnans, we replace with the clean data
                channellabel = cleandata.label{ch};
                channel_tempcleandata = find(strcmp(tempcleandata.label, channellabel));
                cleandata.trial{1,trl}(ch,:) = tempcleandata.trial{1,trl}(channel_tempcleandata,:);
                in=in+1;
            end
        end
    end


    % Finally, we interpolate the bad channels

    load('selected_neighbours.mat'); % Load the neighbours struct

    % NOTE: this struct was created as in the script 'do_group_analysis' and
    % simply copied into the analysis folder

    % we find the electrode positions needed by ft_channelrepair

    elec = ft_read_sens('standard_1020.elc', 'senstype', 'eeg');

    cfg                        = [];
    cfg.method                 = 'nearest';
    cfg.badchannel             = cleandata.label(badchannels);
    cfg.missingchannel         = []; 
    cfg.elec                   = elec;
    cfg.neighbours             = selected_neighbours;
    data_artefact_2            = ft_channelrepair(cfg, cleandata);

    %% Artefact rejection part 3: a final visual artefact rejection

    % We now go through each trial to check if no artefacts remain.
    % Here we can remove some bad trials, but no more bad channels

    cfg                 = [];
    cfg.method          = 'trial';  % Or switch to summary if needed
    cfg.keepchannel     = 'nan';    % when rejecting channels, values are replaced by NaN
    cfg.alim            = [-100, 100];   
    data_cleaned        =  ft_rejectvisual(cfg, data_artefact_2);

    %% Now we find and save rejected trials and channels
    
    badtrial_times2     = data_cleaned.cfg.artfctdef.trial.artifact; % this matrix contains start and end times of each rejected trial
  
    begsample           = [badtrial_times(:,1); badtrial_times2(:,1)];
    endsample           = [badtrial_times(:,2); badtrial_times2(:,2)];
    rejection_round     = [repmat({'Visual rejection 1'}, length(badtrial_times), 1); repmat({'Visual rejection 2'}, length(badtrial_times2), 1)];
    
    badtrials           = table(begsample, endsample, rejection_round);
    
    % For the badchannels, we find the labels corresponding to the bad
    % channels, since it is easier to interpret
    
    badchannels_label   = data_artefact_1.label(badchannels, :);
    
    % Then we save everything
    save(fullfile(output_dir, 'badchannels.mat'), 'badchannels_label');
    save(fullfile(output_dir, 'badchannelnumbers.mat'), 'badchannels');
    save(fullfile(output_dir, 'badtrials.mat'), 'badtrials');
        
    save(fullfile(output_dir, 'data_cleaned.mat'), 'data_cleaned');
        
end % here we end the if loop to do artefact rejection

%% Re-referencing

if exist([output_dir filesep 'data_cleaned.mat'], 'file')
        load([output_dir filesep 'data_cleaned.mat']);
else 
    error('artefact rejected data in data_cleaned cannot be found: either artefact rejection has not been done or it has not been saved properly');
end

cfg                     = [];
cfg.channel             = 'all';  
cfg.reref               = 'yes';
cfg.implicitref         = 'TP9';            % the implicit (non-recorded) reference channel is added to the data representation
cfg.refchannel          = {'TP9', 'TP10'}; % the average of these channels is used as the new reference, note that channel corresponds to the right mastoid (M2)
data_cleaned            = ft_preprocessing(cfg, data_cleaned);

%% Calculate the ERPs for expected (bee) and unexpected (cue) stimuli

% We first perform timelockanalysis on those trials belonging to the
% expected condition (i.e. the bees)

cfg                = [];
cfg.trials         = find(ismember(string(data_cleaned.trialinfo{:,2}), 'bee'));
expected           = ft_timelockanalysis(cfg, data_cleaned);

% And to those of the unexpected condition (i.e. the cues)
cfg = [];
cfg.trials         = find(ismember(string(data_cleaned.trialinfo{:,2}), {'update-cue', 'no-update-cue'}));
unexpected         = ft_timelockanalysis(cfg, data_cleaned);

% And we plot the ERP's
cfg                = [];
cfg.layout         = 'EEG1010.lay';
cfg.interactive    = 'yes';
cfg.showoutline    = 'yes';
cfg.showlabels     = 'yes';
ft_multiplotER(cfg, expected, unexpected)

% Finally we save the data
save(fullfile(output_dir, 'timelock_expected.mat'), 'expected');
save(fullfile(output_dir, 'timelock_unexpected.mat'), 'unexpected');
savefig(gcf, fullfile(output_dir, 'topoplot_expected_unexpected'));

%% Calculate the ERPs for expected (bee) stimuli over number of repetitions

% Initiate vectors
repetition = zeros(size(data_cleaned.trialinfo, 1),1);
count = 0;

for tr = 1:size(data_cleaned.trialinfo, 1)    
    if tr < size(data_cleaned.trialinfo, 1)
        if ismember(string(data_cleaned.trialinfo{tr,2}), 'bee') && ismember(string(data_cleaned.trialinfo{tr + 1,2}), 'bee')
            % this bee is followed by another bee of the same type        
            count = count+1;
            repetition(tr,1) = count; 
        elseif ismember(string(data_cleaned.trialinfo{tr,2}), 'bee') && ~ismember(string(data_cleaned.trialinfo{tr + 1,2}), 'bee')
            % This is the last bee of this repetition of bees
            count = count+1;
            repetition(tr,1) = count;
            count = 0; % we reset count back 
        end
    elseif ismember(string(data_cleaned.trialinfo{tr,2}), 'bee')
        % This is the last bee of the experiment
        count = count+1;
        repetition(tr,1) = count;
    end
end

% Now repetition contains the "repetition number" of the bee, i.e. how many
% times this specific stimulus has been shown uninterrupted by unexpted stimuli 

cfg                = [];
cfg.trials         = find(repetition==1);
repetition1        = ft_timelockanalysis(cfg, data_cleaned); 
%repetition1 contains the average of trials where the bee is shown for the first time

cfg                = [];
cfg.trials         = find(repetition==2);
repetition2        = ft_timelockanalysis(cfg, data_cleaned);
%repetition2 contains the average of trials where the bee is shown for the second time

cfg                = [];
cfg.trials         = find(repetition==3);
repetition3        = ft_timelockanalysis(cfg, data_cleaned);
%repetition3 contains the average of trials where the bee is shown for the third time

% And we plot the ERP's
cfg                = [];
cfg.layout         = 'EEG1010.lay';
cfg.interactive    = 'yes';
cfg.showoutline    = 'yes';
cfg.showlabels     = 'yes';
ft_multiplotER(cfg, repetition1, repetition2, repetition3)

% Now we save the data
save(fullfile(output_dir, 'timelock_repetition1.mat'), 'repetition1');
save(fullfile(output_dir, 'timelock_repetition2.mat'), 'repetition2');
save(fullfile(output_dir, 'timelock_repetition3.mat'), 'repetition3');
savefig(gcf, fullfile(output_dir, 'topoplot_repetitions_expected'))


close all
