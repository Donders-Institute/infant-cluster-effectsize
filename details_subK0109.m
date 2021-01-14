%% Input of the dataset

subjectstr      = mfilename;
sub             = erase(subjectstr, 'details_sub');
sub             = ['sub-' sub];

user = getenv('USER');
if isempty(user)
  user = getenv('UserName');
end

switch user
  case 'Didi'    
    data_dir    = 'C:\Users\Didi\Documents\GitHub\Donders Datasets\BeeG dataset\BIDS';
    output_dir  = ['C:\Users\Didi\Documents\GitHub\Donders Datasets\BeeG dataset\results' filesep sub];    
  case 'roboos'
    data_dir    = '/Volumes/Samsung T3/data/di.dcc.DSC_2020.00134_473/bids';
    output_dir  = ['/Volumes/Samsung T3/data/di.dcc.DSC_2020.00134_473/results' filesep sub];
  otherwise
    error('you have to specify the local directories of the data and this code');
end

addpath(data_dir)
cd(data_dir)

if ~exist(output_dir, 'dir')
   mkdir(output_dir);    
end

%% Find badchannels and badtrials

% If artefact rejection has been done, we can extract the rejected channels
% and trials, so they can easily be read here

% First the trials, they should be in the participants.tsv

participants                       = read_tsv([data_dir filesep 'participants.tsv']);

if ~isnan(participants.excluded_trials(subjectnumber, 1))
    % We filled out a number here: analysis has been done
    number_trials_rejected         = participants.excluded_trials(subjectnumber, 1);
    % Then the trials that have been rejected are stored in a mat struct in
    % the output dir folder
    if exist([output_dir filesep 'badtrials.mat'], 'file')==2
        trialrejection             = load([output_dir filesep 'badtrials.mat']); 
    end
    
    % We don't have to do the artefact rejection again
    do_artefact_rejection_test1    = 0;
end

% Now for the bad channels, they are in the channels tsv
channels_folder                    = dir([data_dir filesep sub filesep 'eeg' filesep '*channels*.tsv']);
channels_name                      = [data_dir filesep sub filesep 'eeg' filesep channels_folder.name];
channels_tsv                       = read_tsv(channels_name);

if sum(ismember(channels_tsv.Properties.VariableNames, 'status'))~=0
    % the channels have been analyzed
    if exist([output_dir filesep 'badchannels.mat'], 'file')==2
        channelrejection            = load([output_dir filesep 'badchannels.mat']); 
    end
    % We don't have to do artefact rejection again
    do_artefact_rejection_test2    = 0;
end

%% We test if we want to do artefact rejection again

if do_artefact_rejection_test1 == 0 && do_artefact_rejection_test2 == 0
    disp('')
    disp('The following trials have been rejected')
    printstruct(trialrejection.badtrials)
    disp('')
    disp('The following channels have been rejected')
    printstruct(channelrejection.badchannels_label)
    yn=input('Press y to continue and n to redo artefact rejection','s');
    if strcmp(yn,'y')==1
        do_artefact_rejection = 0;
    elseif strcmp(yn,'n')==1
        do_artefact_rejection = 1;
    end
end
