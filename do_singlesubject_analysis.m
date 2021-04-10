%% Analysis script to analyze a single subject

% Set paths
do_setpath

% Display current subject being processed
fprintf('\n')
disp('------------------------------------')
disp (['Doing analysis for subject: ' sub])
disp('------------------------------------')
fprintf('\n')

% Specifying results directory for this specific subject
output_dir = fullfile(results, sub);

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%%  Test whether artefact rejection has already been performed or not

% Check if the badtrials and badchannels .mat files already exist
if exist([output_dir filesep 'badtrials.mat'], 'file') && exist([output_dir filesep 'badchannels.mat'], 'file')
    % The analysis has alreay been performed: ask the user if it needs
    % to be repeated
    yn = input('Artefact rejection has already been done. Do you want to re-do it? [press y / n]','s');
    if strcmp(yn,'y')
        do_artefact_rejection = 1;
    elseif strcmp(yn,'n')
        do_artefact_rejection = 0;
    end
else
    % If artefact rejection has not been done yet: specify that it will be
    % performed
    do_artefact_rejection = 1;
end


%% Find the data and headerfiles

cfg         = [];
cfg.dataset = [bidsroot filesep sub filesep 'eeg' filesep sub '_task-audiovisual_eeg.vhdr' ];
hdr         = ft_read_header(cfg.dataset);

%% Now extract the trials from the events.tsv

t              = readtable([bidsroot filesep sub filesep 'eeg' filesep sub '_task-audiovisual_events.tsv'], 'FileType', 'text');
trl            = t( : , 3:8);

% The variable raw now contains the eeg tsv info:
% Column 1: begsample
% Column 2: endsample
% Column 3: offset
% Column 4: marker number
% Column 5: stimulus
% Column 6: location of the stimulus on the screen

% Then, re-define trials. For the purpose of this, we focus on standard
% (bees) and oddball (cues) stimuli,
% and take a time window of -500 ms and +1000 ms with respect to stimulus
% onset

stimuli             = {'bee', 'update-cue', 'no-update-cue'};
str                 = string(trl{:,5});
index_stimuli       = ismember(str, stimuli);
trl_cueonly         = trl(index_stimuli, :); % now only standeard (bee) and oddball (cue) trials are left

pre_stim_samples    = round(0.5 * hdr.Fs); % The pre-stim period is 0.5 s
post_stim_samples   = round(1 * hdr.Fs);   % The post-stim period is 1 s

% If the first stimulus does not have enough pre-stimulus time for the required window, remove it
if trl_cueonly{1,1} < pre_stim_samples
    trl_cueonly     = trl_cueonly(2:end, :);
    trl_new         = trl_cueonly;
else
    trl_new         = trl_cueonly;
end

trl_new{:,1}        = trl_cueonly{:,1}-pre_stim_samples; % subtract the pre-stimulus samples from begsample (onset stimulus) to find the start of the time window of interest
trl_new{:,2}        = trl_cueonly{:,1}+post_stim_samples; % add the post stim samples to begsample (onset stimulus) to find the end of the time window of interest
trl_new{:,3}        = - pre_stim_samples;

% Then add the new trials to cfg
cfg.trl             = trl_new;

% Save the new trials
save(fullfile(output_dir, 'trials.mat'), 'trl_new');

%% Pre-processing on the epochs that are defined as trials

% We perform pre-processing similarly to the published work on this study:
% Kayhan et al. Developmental Cognitive Neuroscience, 2019
% That is: 5 sec padding for high pass filtering at 1 Hz and baseline
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
% 1) A first pass visual rejection of large artifacts
% 2) ICA to detect and correct for artefacts like those caused by
% eye-movements
% 3) A second pass of visual artefact rejection to remove any remaining artefacts
% 4) Interpolating channels
% This is only necessary once, if it has been done, we skip this step


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define neighbours: 
% To interpolate neighbouring channels and to run
% cluster-based statistic including channels first find neighbouring
% channels for the current recoding layout.
% This was done once with the code below. Note, that this is now commented out as this struct only needs to created once
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% To start read in the channels used in the EEG study

% label = expected.cfg.channel;
%
% % Read sensor positions from the standard 1020 3D format available on fieldtrip
% elec = ft_read_sens('standard_1020.elc', 'senstype', 'eeg');
%
% if exist('newlayout')
%     clear newlayout
% end
%
% % Now loop through all used channels, to create a new layout struct
% % containing only those channels used in this particular study
%
% for ii = 1:length(label)
%       index                     = strcmp(label{ii}, elec.label);
%       newlayout.chanpos(ii, :)  = elec.chanpos(index, :);
%       newlayout.chantype{ii, 1} = elec.chantype{index};
%       newlayout.chanunit{ii, 1} = elec.chanunit{index};
%       newlayout.elecpos(ii, :)  = elec.elecpos(index, :);
%       newlayout.label{ii, 1}    = elec.label{index};
%       newlayout.type            = elec.type;
%       newlayout.unit            = elec.unit;
% end
%
% % Then plot this in 3D
%
% fig = ft_plot_sens(newlayout, 'label', 'label');
%
% %% Draw lines in the figure to connect channels and create a neighbours struct
%
% % Loop through all channels
%
% for cc = 1:length(newlayout.label)
%
%      % Find the corresponding label
%      channel_to_connect = newlayout.label{cc};
%
%      % Allow the user to provide a cell array of channels to connect
%      % the selected channel to
%
%      fprintf('\n')
%      disp(['channel to connect = ' channel_to_connect ', channel number ' num2str(cc) ' of ' num2str(length(newlayout.label))]);
%      text = 'provide a row cell array named channels to connect to channel to connect: ';
%      channels = input(text);
%
%      % Create a line between the selected channel and the channels
%      % to connect that channel to
%
%      connect    = find(strcmp(newlayout.label, channel_to_connect));
%
%      for ci = 1:length(channels)
%          connect1    = find(strcmp(newlayout.label, channels(ci)));
%          linepoints  = [newlayout.chanpos(connect, :); newlayout.chanpos(connect1, :)];
%          line(linepoints(:,1), linepoints(:,2), linepoints(:,3));
%      end
%
%      % Create the neighbours struct
%
%      selected_neighbours(cc).label       = channel_to_connect;
%      selected_neighbours(cc).neighblabel = channels;
%
%  end
%
% % Save the neighbours struct and created image
% save(fullfile(scripts, 'selected_neighbours.mat'), 'selected_neighbours');
% savefig(gcf, fullfile(scripts, 'selected_neighbours_plot'));
%
% % When this step of the code is finished, comment it out to continue with
% % manual optimization in the next step
%
% %%  Last step: make some manual changes if necessary
%
% if exist('selected_neighbours.mat')
%     load('selected_neighbours.mat');
% end
%
% if exist('selected_neighbours_plot.fig')
%     openfig('selected_neighbours_plot.fig');
% end
%
% % Add a line between two points that were missed by the previous step:
%
% connect     = find(strcmp(newlayout.label, 'Fp1'));
% connect1    = find(strcmp(newlayout.label, 'Fp2'));
% linepoints  = [newlayout.chanpos(connect, :); newlayout.chanpos(connect1, :)];
% line(linepoints(:,1), linepoints(:,2), linepoints(:,3));
%
% % And add the connection to the neighours struct
% selected_neighbours(1).neighblabel  = [selected_neighbours(1).neighblabel, {'Fp2'}];
% selected_neighbours(2).neighblabel  = [selected_neighbours(2).neighblabel, {'Fp1'}];
%
% % Then do a sanity check, make sure that neighbours are fully bidirectional
% % I.e. if Cz has FCz as a neighbour, then FCz should have Cz as a neighbour
%
% for ii = 1:length(selected_neighbours)
%     channel = selected_neighbours(ii).label;
%     for tt  = 1:length(selected_neighbours(ii).neighblabel)
%         idx = find(strcmp(label, selected_neighbours(ii).neighblabel{tt}));
%         if sum(strcmp(selected_neighbours(idx).neighblabel, channel))==0
%           selected_neighbours(idx).neighblabel = [selected_neighbours(idx).neighblabel, {channel}];
%         end
%     end
% end
%
% % And replot so that all lines are now corrected
%
% close all
%
% fig = ft_plot_sens(newlayout, 'label', 'label');
%
% for ii = 1:length(selected_neighbours)
%     channel         = selected_neighbours(ii).label;
%     connect         = find(strcmp(newlayout.label, channel));
%     for tt          = 1:length(selected_neighbours(ii).neighblabel)
%         connect1    = find(strcmp(newlayout.label, selected_neighbours(ii).neighblabel{tt}));
%         linepoints  = [newlayout.chanpos(connect, :); newlayout.chanpos(connect1, :)];
%         line(linepoints(:,1), linepoints(:,2), linepoints(:,3));
%     end
% end
%
% % Finally, save this plot and new neighbours struct
%
% save(fullfile(scripts, 'selected_neighbours.mat'), 'selected_neighbours');
% savefig(gcf, fullfile(scripts, 'selected_neighbours_plot'));


if do_artefact_rejection==1 % if the artefact rejection has not been conducted yet
    
    %% Artefact rejection part 1: A first pass visual rejection using ft_rejectvisual
    
    % Start with trial (or summary) view to get an overview of the data, and to reject channels or trials that show large artefacts
    
    cfg                 = [];
    cfg.method          = 'trial';  % Or switch to 'summary' if needed
    cfg.keepchannel     = 'nan';    % when rejecting channels, values are replaced by NaN
    cfg.ylim            = [-100, 100];
    data_artefact_1     =  ft_rejectvisual(cfg, data);
       
    % Then find the rejected trials
    badtrial_times      = data_artefact_1.cfg.artfctdef.trial.artifact; % this matrix contains start and end times of each rejected trial
    
    %% Artefact rejection part 2: ICA
    
    % Before doing ICA remove bad channels
    badchannels             = find(all(isnan(data_artefact_1.trial{1}), 2));
    channels                 = data_artefact_1.label;
    channels(badchannels, :) = [];
    
    % Then set up the cfg struct and perform ICA
    
    cfg                      = [];
    cfg.channel              = channels;
    cfg.method               = 'runica'; 
    comp                     = ft_componentanalysis(cfg, data_artefact_1);
    
    % Then  plot the first 20 components
    cfg                      = [];
    cfg.component            = 1:20; % We plot the first 20 components
    cfg.layout               = 'EEG1010.lay';
    cfg.comment              = 'no';
    cfg.zlim                 = 'maxabs';
    
    ft_topoplotIC(cfg, comp)
    
    disp('')
    input('Press Enter to continue')
    close all
    
    % Plot the time course of all components again to check whether the components
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
    
    % Also exclude all bad channels specified during the first visual artefact rejection (containing NaNs)
    exclude                      = find(~ismember(data_artefact_1.label, channels));
    tempdata                     = data_artefact_1;
    for trl=1:size(tempdata.trial, 2)
        tempdata.trial{1,trl}(exclude,:) = [];
        tempdata.label = channels;
    end
    
    cfg                          = [];
    cfg.component                = rejcom; % to be removed components
    cfg.channel                  = channels;
    tempcleandata                = ft_rejectcomponent(cfg, comp, tempdata);
    
    % Then copy everything back into the full data structure with the 'bad
    % channels' as nans
    
    ica_cleandata                    = data_artefact_1;
    for trl = 1:size(ica_cleandata.trial, 2) % Loop through trials
        in=1;
        for ch = 1:size(ica_cleandata.trial{1,trl},1) % Loop through channels
            if ~isnan(ica_cleandata.trial{1,trl}(ch,1)) % if the channel was not excluded previously replace the original data with the ICA-cleaned data
                channellabel = ica_cleandata.label{ch};
                channel_tempcleandata = find(strcmp(tempcleandata.label, channellabel));
                ica_cleandata.trial{1,trl}(ch,:) = tempcleandata.trial{1,trl}(channel_tempcleandata,:);
                in=in+1;
            end
        end
    end
 
        
    %% Artefact rejection part 3: A second pass of visual artefact rejection
    
    % Go through each trial to check if any artefacts remain and remove them
    
    cfg                 = [];
    cfg.method          = 'trial';  % Or switch to summary if needed
    cfg.keepchannel     = 'nan';    % when rejecting channels, values are replaced by NaN
    cfg.ylim            = [-100, 100];
    data_cleaned        =  ft_rejectvisual(cfg, ica_cleandata);
    
    %% Interpolate bad channel and save info on rejected trials and channels
    
    badtrial_times2     = data_cleaned.cfg.artfctdef.trial.artifact; % this matrix contains start and end times of each rejected trial
    
    begsample           = [badtrial_times(:,1); badtrial_times2(:,1)];
    endsample           = [badtrial_times(:,2); badtrial_times2(:,2)];
    rejection_round     = [repmat({'Visual rejection 1'}, size(badtrial_times, 1), 1); repmat({'Visual rejection 2'}, size(badtrial_times2,1), 1)];
    
    badtrials           = table(begsample, endsample, rejection_round);
          
    % Find the bad channels
    badchannels         = find(all(isnan(data_cleaned.trial{1}), 2));
    
    % Find the corresponding label for each bad channel 
    badchannels_label   = data_cleaned.label(badchannels, :);
    
    load(fullfile(scripts, 'selected_neighbours.mat')); % Load the neighbours struct
       
    % Then find the electrode positions needed by ft_channelrepair
    
    elec = ft_read_sens('standard_1020.elc', 'senstype', 'eeg');
    
    cfg                        = [];
    cfg.method                 = 'weighted';
    cfg.badchannel             = data_cleaned.label(badchannels);
    cfg.missingchannel         = [];
    cfg.elec                   = elec;
    cfg.neighbours             = selected_neighbours;
    data_cleaned               = ft_channelrepair(cfg, data_cleaned);
    
    % Then save everything
    save(fullfile(output_dir, 'badchannels.mat'), 'badchannels_label');
    save(fullfile(output_dir, 'badchannelnumbers.mat'), 'badchannels');
    save(fullfile(output_dir, 'badtrials.mat'), 'badtrials');
    
    save(fullfile(output_dir, 'data_cleaned.mat'), 'data_cleaned');
    
end % end of artefact rejection

%% Re-referencing

if exist([output_dir filesep 'data_cleaned.mat'], 'file')
    load([output_dir filesep 'data_cleaned.mat']);
else
    error('Artefact rejected data in data_cleaned cannot be found: either artefact rejection has not been done or it has not been saved properly');
end

cfg                     = [];
cfg.channel             = 'all';
cfg.reref               = 'yes';
cfg.implicitref         = 'TP9';            % the implicit reference channel is added to the data representation
cfg.refchannel          = {'TP9', 'TP10'}; % the average of these channels is used for re-reference
data_cleaned            = ft_preprocessing(cfg, data_cleaned);

%% Calculate ERPs for standard (bee) and oddball (cue) stimuli
no_trls_left = [];

% Perform timelockanalysis on standard (bee) stimuli
cfg                = [];
cfg.trials         = find(ismember(string(data_cleaned.trialinfo{:,2}), 'bee'));
if ~isempty(cfg.trials) % if there are trials in this condition
    standard           = ft_timelockanalysis(cfg, data_cleaned);
    % Save the data
    save(fullfile(output_dir, 'timelock_standard.mat'), 'standard');
else
    no_trls_left = 1;
end

% Perform timelockanalysis on oddball (cue) stimuli
cfg = [];
cfg.trials          = find(ismember(string(data_cleaned.trialinfo{:,2}), {'update-cue', 'no-update-cue'}));
if ~isempty(cfg.trials) % if there are trials in this condition
    oddball             = ft_timelockanalysis(cfg, data_cleaned);
    % Save the data
    save(fullfile(output_dir, 'timelock_oddball.mat'), 'oddball');
else
    no_trls_left = 1;
end

% Plot ERP's if there are data for both conditions
if isempty(no_trls_left)
    cfg                = [];
    cfg.layout         = 'EEG1010.lay';
    cfg.interactive    = 'yes';
    cfg.showoutline    = 'yes';
    cfg.showlabels     = 'yes';
    ft_multiplotER(cfg, standard, oddball)
    
    % Save the figure
    savefig(gcf, fullfile(output_dir, 'topoplot_standard_oddball'));
end
close all
