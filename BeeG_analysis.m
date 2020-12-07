%% Analysis script for the BeeG dataset

%% Input of the dataset

clear;

user = getenv('USER');
if isempty(user)
  user = getenv('UserName');
end

switch user
  case 'Didi'    
    bidsroot    = 'C:\Users\Didi\Documents\GitHub\Donders Datasets\BeeG dataset\BIDS';
  case 'roboos'
    bidsroot    = '/Volumes/Samsung T3/data/di.dcc.DSC_2020.00134_473/bids';
  otherwise
    errror('you have to specify the local directories of the data and this code');
end

addpath(bidsroot)
cd(bidsroot)

%% Find the data and headerfiles

% Read the participants tsv to find subject info
[data, header, raw] = tsvread([bidsroot filesep 'participants.tsv']);
sub = raw(2:end, 1);

% Then we start looping over subjects
ii = 1;

% And find the headerfile and eeg file

cfg.dataset = [bidsroot filesep sub{ii} filesep 'eeg' filesep sub{ii} '_task-audiovisual_eeg.vhdr' ];
hdr         = ft_read_header(cfg.dataset);


%% Now we can extract the trials from the events.tsv

[data, header, raw] = tsvread([bidsroot filesep sub{ii} filesep 'eeg' filesep sub{ii} '_task-audiovisual_events.tsv']);

% We don't need begsample or endsample
trl                 = raw(2:end, 3:8);


% Then re-define trials. We are only interested in bees and cues
% And take a time window of -500 ms and + 1000 ms

stimuli             = {'bee', 'update-cue', 'no-update-cue'};
str                 = string(trl(:,5));
index_stimuli       = ismember(str, stimuli);
trl_cueonly         = trl(index_stimuli, :); % now only bee and cue trials are left

pre_stim_samples    = round(0.5 * hdr.Fs); % The pre-stim period is 0.5 s
post_stim_samples   = round(1 * hdr.Fs);   % The post-stim period is 1 s
trl_new             = trl_cueonly;
trl_new(:,1)        = num2cell(str2double(trl_cueonly(:,1))-pre_stim_samples); % We extract the pre stim samples from begsample (the first colum of trl_cueonly) to find the new begsample
trl_new(:,2)        = num2cell(str2double(trl_cueonly(:,1))+post_stim_samples); % We add the post stim samples to begsample to find the new endsample
trl_new(:,3)        = num2cell(str2double(trl_cueonly(:,3)));

% Then add the new trials to cfg
cfg.trl             = trl_new;

%% Now, we perform pre-processing only on the episodes that are defined as trials

% Baseline-correction options
cfg.demean          = 'yes';
cfg.baselinewindow  = [-0.2 0];

% Fitering options
cfg.lpfilter        = 'yes';
cfg.lpfreq          = 100;

% Re-referencing options - see explanation above
cfg.implicitref   = 'LM';
cfg.reref         = 'yes';
cfg.refchannel    = {'LM' 'RM'};

data = ft_preprocessing(cfg);

%% Now let's look at the data with databrowser

cfg = [];  % use only default options
ft_databrowser(cfg, data);

%% Compare this to the unprocessed data

cfg         = [];
cfg.dataset = 's04.vhdr';
cfg.viewmode = 'vertical';
ft_databrowser(cfg);

%% Next: extracting the EOG signal

% EOGV channel
cfg              = [];
cfg.channel      = {'53' 'LEOG'}; % for the vertical EOG
cfg.reref        = 'yes';
cfg.implicitref  = []; % this is the default, we mention it here to be explicit
cfg.refchannel   = {'53'};
eogv             = ft_preprocessing(cfg, data); % eogv now contains a referenced LEOG channel (LEOG-53)


% only keep one channel, and rename to eogv
cfg              = [];
cfg.channel      = 'LEOG';
eogv             = ft_selectdata(cfg, eogv);
eogv.label       = {'eogv'};

% EOGH channel, the same for the horizontal referencing
cfg              = [];
cfg.channel      = {'57' '25'};
cfg.reref        = 'yes';
cfg.implicitref  = []; % this is the default, we mention it here to be explicit
cfg.refchannel   = {'57'};
eogh             = ft_preprocessing(cfg, data);

% only keep one channel, and rename to eogh
cfg              = [];
cfg.channel      = '25';
eogh             = ft_selectdata(cfg, eogh);
eogh.label       = {'eogh'};

% only keep all non-EOG channels
cfg         = [];
cfg.channel = setdiff(1:60, [53, 57, 25]);        % you can use either strings or numbers as selection
data        = ft_selectdata(cfg, data);

% append the EOGH and EOGV channel to the 60 selected EEG channels
cfg  = [];
data = ft_appenddata(cfg, data, eogv, eogh);

%% Check the data again

cfg = [];  % use only default options
ft_databrowser(cfg, data);

%% Let's set up the layout of the channels

cfg        = [];
cfg.layout = 'mpi_customized_acticap64.mat'; % they in this .mat file, fieldtrip has a number of default formats
ft_layoutplot(cfg);

%% Then artifact rejection

% Here we do a visual rejection, not automated. First in channel mode.

cfg        = [];
cfg.method = 'channel';
ft_rejectvisual(cfg, data) % here, we do not assign a new value to is, so rejects will not be saved.

% To save, do data_artifact_rej = ft_rejectvisual(cfg, data)

%% Artifact rejection in summary mode, still visual

cfg          = [];
cfg.method   = 'summary';
cfg.layout   = 'mpi_customized_acticap64.mat';  % for plotting individual trials
cfg.channel  = [1:60];                          % do not show EOG channels
data_clean   = ft_rejectvisual(cfg, data);

%% Now check the cleaned data

cfg = [];
cfg.viewmode = 'vertical';
ft_databrowser(cfg, data_clean);

%% Calculate the ERPs

% use ft_timelockanalysis to compute the ERPs
cfg = [];
cfg.trials = find(data_clean.trialinfo==1); % these are the trials of one condition (1)
task1 = ft_timelockanalysis(cfg, data_clean); % ft_timelockanalysis computes the timelocked average ERP/ERF

cfg = [];
cfg.trials = find(data_clean.trialinfo==2); % these are the trials of the other condition (2)
task2 = ft_timelockanalysis(cfg, data_clean); 

cfg = [];
cfg.layout = 'mpi_customized_acticap64.mat';
cfg.interactive = 'yes';
cfg.showoutline = 'yes';
ft_multiplotER(cfg, task1, task2)

%% Then look at the difference ERP waves between task1 and task2

cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
difference = ft_math(cfg, task1, task2);

cfg = [];
cfg.layout      = 'mpi_customized_acticap64.mat';
cfg.interactive = 'yes';
cfg.showoutline = 'yes';
ft_multiplotER(cfg, difference);


