%% Conversion of the original "sourcedata" into BIDS.
%
% These scripts and the data in BIDS format are part of Meyer, M., Lamers, D., Kayhan,
% E., Hunnius, S., & Oostenveld, R. (2021) Fostering reproducibility in developmental
% EEG research by using BIDS, cluster-based permutation tests and reporting
% effectsizes (in preparation)
%
% The infant EEG dataset is originally described in Kayhan, E., Meyer, M., O'Reilly,
% J. X., Hunnius, S., & Bekkering, H. (2019). Nine-month-old infants update their
% predictive models of a changing environment. Developmental cognitive neuroscience,
% 38, 100680.)
%
% The code in this script is referred to as script section 1

%% 1.1 Specification of folders

clear

do_setpath

cd(sourcedata)

% Delete the current BIDS folder if it already exists to have create the latest version
if exist(bidsroot, 'dir')
  rmdir(bidsroot, 's');
end

%% 1.2 Subject information

% Read the excel file containing subject info
subject_file            = [sourcedata filesep 'S_BeeG_participant_info_all'];
[age, text]             = xlsread(subject_file);
sub                     = text(2:end, 1);
sex                     = text(2:end, 3);

%% 1.3 General information for the data2bids function

% Start looping over subjects
for ii = 1:length(sub)
  
  % The data is in BrainVision format, so the file does not need to be
  % adjusted and can simply be copied
  cfg                                         = [];
  cfg.method                                  = 'copy';
  cfg.bidsroot                                = bidsroot;
  cfg.datatype                                = 'eeg';
  cfg.writejson                               = 'replace';
  
  %% 1.4 Create the metafile in the dataset: dataset_description.json
  
  cfg.dataset_description.Name                = 'Nine-month-old infants update their predictive models of a changing environment';
  cfg.dataset_description.DatasetType         = 'raw';
  cfg.dataset_description.BIDSVersion         = '1.2.0';
  cfg.dataset_description.Authors             = {'Ezgi Kayhan','Marlene Meyer', 'Jill X O''Reilly', 'Sabine Hunnius', 'Harold Bekkering'};
  cfg.dataset_description.License             = 'ODC-ODbL-1.0';
  
  % Other options could be the following
  %   cfg.dataset_description.Acknowledgements    = string
  %   cfg.dataset_description.HowToAcknowledge    = string
  %   cfg.dataset_description.Funding             = string or cell-array of strings,
  %   cfg.dataset_description.ReferencesAndLinks  = string or cell-array of strings
  %   cfg.dataset_description.DatasetDOI          = string
  
  %% 1.5 Create a file with project level participant information: participants.tsv
  
  cfg.sub                               = sub{ii};
  cfg.sub(cfg.sub=='_')                 = []; % remove underscores from the subject identifier
  
  cfg.participants.age                  = age(ii);
  cfg.participants.sex                  = sex{ii};
  
  %% 1.6 Select the indiviudal EEG dataset per participant
  
  cfg.dataset                           = ['Raw_data_infant' filesep  sub{ii} '.vhdr'];
  if cfg.dataset
    hdr                                   = ft_read_header(cfg.dataset);
  end
  
  %% 1.7 Create metadata information on the EEG recordings
  
  % Describing the task
  cfg.TaskName                                       = 'audiovisual';
  cfg.TaskDescription                                = 'Infants observed a sequence of repeated audio-visual stimuli, followed by a cue signalling subsequent updating or no-updating of the sequence. Respectively, the sequence continued with the same (no-update)or new (updated) stimuli.';
  cfg.Instructions                                   = 'Parents were instructed to keep the interaction with their infant minimal during the measurement, infants received no instructions';
  
  % Describing the recording setup
  cfg.InstitutionName                                = 'The Donders Institute for Brain, Cognition and Behaviour';
  cfg.InstitutionAddress                             = 'Heyendaalseweg 135, 6525 AJ Nijmegen, the Netherlands';
  cfg.InstitutionalDepartmentName                    = 'Donders Centre for Cognition';
  
  cfg.Manufacturer                                   = 'Brain Products GmbH';
  cfg.ManufacturersModelName                         = 'BrainAmp Standard';
  cfg.eeg.CapManufacturer                            = 'Brain Products GmbH';
  cfg.eeg.CapManufacturersModelName                  = 'actiCAP 32Ch';
  cfg.eeg.EEGPlacementScheme                         = '10-20';
  cfg.eeg.EEGReference                               = 'TP9';
  cfg.eeg.EEGGround                                  = 'AFz';
  cfg.eeg.SamplingFrequency                          = 500;
  
  % NOTE: the amplifier always samples at 5000 Hz in hardware, the data
  % is then downsampled to 500 Hz in software
  
  cfg.eeg.PowerLineFrequency                         = 50;
  cfg.eeg.HardwareFilters.LowCutoff.Frequency        = 0.1;
  cfg.eeg.HardwareFilters.HighCutoff.Frequency       = 1000;
  cfg.eeg.SoftwareFilters.LowCutoff.Frequency        = 0.1;
  cfg.eeg.SoftwareFilters.HighCutoff.Frequency       = 200;
  cfg.eeg.EEGChannelCount                            = 32;
  
  % Describing the recording
  cfg.eeg.RecordingType                              = 'continuous';
  % cfg.eeg.RecordingDuration                        = % Read automatically
  % cfg.eeg.EpochLength                              = % Read automatically
  
  %% 1.8 Create metadata on markers/triggers used in the EEG recordings: events.tsv
  
  % To do this, create segment data into events using ft_define_trial
  cfg_trials                      = cfg;
  cfg_trials.trialdef.eventtype   = 'Stimulus';
  trl                             = trialfun_BeeG(cfg_trials);
  cfg.events                      = trl;
  
  %% 1.9 Create metadata on the channels used in the EEG recordings: channels.tsv
  
  % Double info with eeg.tsv --> here only fill it out if it is channel specific
  cfg.channels.name               = hdr.label;
  cfg.channels.type               = repmat({'EEG'}, 32, 1);  % Type of channel
  cfg.channels.units              = repmat({'uV'}, 32, 1);% Physical unit of the data values recorded by this channel in SI
  cfg.channels.sampling_frequency = repmat(500, 32, 1); % Sampling rate of the channel in Hz.
  cfg.channels.low_cutoff         = repmat(0.1, 32, 1); % Frequencies used for the hardware high-pass filter applied to the channel in Hz
  cfg.channels.high_cutoff        = repmat(1000, 32, 1); % Frequencies used for the hardware low-pass filter applied to the channel in Hz.
  cfg.channels.notch              = repmat(nan, 32, 1); % Frequencies used for the notch filter applied to the channel, in Hz. If no notch filter applied, use n/a.
  
  % Other options could be the following
  % cfg.channels.software_filters   = {' "Low Cutoff": 0.1', ' "High Cutoff": 125'}; % List of temporal and/or spatial software filters applied.
  % cfg.channels.description        = % OPTIONAL. Brief free-text description of the channel, or other information of interest. See examples below.
  % cfg.channels.status             = % OPTIONAL. Data quality observed on the channel (good/bad). A channel is considered bad if its data quality is compromised by excessive noise. Description of noise type SHOULD be provided in [status_description].
  % cfg.channels.status_description = % OPTIONAL. Freeform text description of noise or artifact affecting data quality on the channel. It is meant to explain why the channel was declared bad in [status].
  
  
  %% Call data2bids
  
  data2bids(cfg);
  
  %% 1.10 Create a description of the markers/triggers used in the EEG recordings: events.json
  
  events_json                                 = [];
  events_json.onset.description               = 'Onset of the event';
  events_json.onset.units                     = 'seconds';
  events_json.duration.description            = 'Duration of the event';
  events_json.duration.units                  = 'seconds';
  events_json.begsample.description           = 'Sample at which event begins (measured from start of recording)';
  events_json.begsample.units                 = 'sample number';
  events_json.endsample.description           = 'Sample at whcih event ends (measured from start of recording)';
  events_json.endsample.units                 = 'sample number';
  events_json.offset.description              = 'Offset from begsample till start of the trial';
  events_json.offset.units                    = 'sample number';
  events_json.marker.description              = 'Marker number corresponding to this event as indicated in the .vmrk file';
  events_json.stimulus.description            = 'Type of stimulus presented to the infant';
  events_json.stimulus.levels                 = {'blank screen: presentation of blank screen before fixation cross',...
    'fixation cross: presentation of a fixation cross (see subfoler stimuli\fixation)and a sound(see subfoler stimuli\Ding_Sound_Effect)',...
    'bee: presentation of a bee image in the same location and of the same colour as the previous bee (see subfoler stimuli\bee) along with a jumping sound (see subfolder stimuli\sound1-9)',...
    'post update-cue bee: presentation of a bee image after an update cue, in a different location and of a different colour as the previous bee (see subfoler stimuli\bee) along with a jumping sound (see subfolder stimuli\sound 1-9)',...
    'post no-update-cue bee: presentation of a bee image after a no-update cue, in the same location and of the same colour as the previous bee (see subfoler stimuli\bee) along with a jumping sound (see subfolder stimuli\sound 1-9)',...
    'update-cue: presentation of the update-cue (see subfolder stimuli\circle or triangle) and a sound (see subfolder stimuli\sound10 or sound11)',...
    'no-update-cue: presentation of the no-update-cue (see subfolder stimuli\circle or triangle) and a sound (see subfolder stimuli\sound10 or sound11)',};
  events_json.location_bee.description         = 'Location of the bee image';
  events_json.location_bee.units               = 'degrees';
  
  foldername                                  = [bidsroot filesep 'sub-' cfg.sub filesep 'eeg'];
  filename                                    = [foldername filesep 'sub-' cfg.sub '_task-' cfg.TaskName '_events.json'];
  
  ft_write_json(filename, events_json);
  
  
end % End subject loop

%% 1.11 Create description of metadata on participants: participants.json

participants_json.participant_id.description    = 'Subject identifier';
participants_json.age.description               = 'age of each subject';
participants_json.age.units                     = 'days';
participants_json.sex.description               = 'gender of each subject';
participants_json.sex.levels                    = {'girl: female', 'boy: male'};

filename                                        = [bidsroot filesep 'participants.json'];

ft_write_json(filename, participants_json);

%% Add the matlab code used to generate BIDS to a subfolder

destination                                     = [bidsroot filesep 'code'];
this_script                                     = [mfilename('fullpath') '.m'];
trialfun_script                                 = [fileparts(which(mfilename)) filesep 'trialfun_BeeG.m'];

mkdir(destination);
copyfile(this_script, destination);
copyfile(trialfun_script, destination);

%% Create a sourcedata folder for the logfiles and a Stimuli folder for the video and audio stimuli

% Let's create a folder for the logfiles
str                                             = [bidsroot filesep 'sourcedata' filesep 'Logfiles'];
mkdir(str)

% And for the video files
strvideo                                        = [bidsroot filesep 'stimuli'];
mkdir(strvideo)

% Then copy the data there
logfiles                                        = [sourcedata filesep 'Logfiles'];
logfiles                                        = dir(logfiles);
videofile                                       = [sourcedata filesep 'Stimuli'];
videofile                                       = dir(videofile);

copyfile([sourcedata filesep 'Logfiles' filesep logfiles(1).name], str);
copyfile([sourcedata filesep 'Stimuli' filesep logfiles(1).name], strvideo);

%% Exlude scans.tsv from bidsvalidator, it does not support the vhdr name of the dataset as valid

destination = fullfile(bidsroot, '.bidsignore');
fileID = fopen(destination,'w');
fprintf(fileID,'*_scans.tsv\n');
fclose(fileID);

%% Make a README and CHANGES file as placeholders

destination = fullfile(bidsroot, 'CHANGES');
fileID = fopen(destination,'w');
fprintf(fileID,'Revision history\n\n\n');
fclose(fileID);

destination = fullfile(bidsroot, 'README');
fileID = fopen(destination,'w');
fprintf(fileID,'\n\n\n');
fclose(fileID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUBFUNCTION

function val = output_compatible(val)
if istable(val)
  fn = val.Properties.VariableNames;
  for i=1:numel(fn)
    val.(fn{i}) = output_compatible(val.(fn{i}));
  end
elseif iscell(val)
  % use recursion to make all elements compatible
  val = cellfun(@output_compatible, val, 'UniformOutput', false);
elseif isnumeric(val) && numel(val)>1 && any(isnan(val))
  % convert and use recursion to make all elements compatible
  val = num2cell(val);
  val = cellfun(@output_compatible, val, 'UniformOutput', false);
else
  % write [] as 'n/a'
  % write nan as 'n/a'
  % write boolean as 'True' or 'False'
  if isempty(val)
    val = 'n/a';
  elseif isnan(val)
    val = 'n/a';
  elseif islogical(val)
    if val
      val = 'True';
    else
      val = 'False';
    end
  end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUBFUNCTION

function tsv_obs = read_excel_observation(excel_observation)

CodedPeriod                                    = cell(length(excel_observation), 1);
for ll = 1:length(excel_observation)
  if excel_observation(ll, 1) == 95
    CodedPeriod(ll, :) = {'fixation cross'};
  elseif excel_observation(ll, 1)>100 && excel_observation(ll, 1)<200
    CodedPeriod(ll, :) = {'intro video'};
  elseif excel_observation(ll, 1)>2000000 && excel_observation(ll, 1)<22000000
    % there are nine phases of this:
    phase = num2str(excel_observation(ll, 1));
    phase_str = ['stimulus video - phase' phase(end)];
    CodedPeriod(ll, :) = {phase_str};
  else
    phase = num2str(excel_observation(ll, 1));
    phase_str = ['peek-a-boo video - phase' phase(end)];
    CodedPeriod(ll, :) = {phase_str};
  end
end
annotation                                    = excel_observation(:, 2);
onset                                         = nan(size(annotation));
duration                                      = nan(size(annotation));
tsv_obs                                       = table(onset, duration, CodedPeriod, annotation);

end
