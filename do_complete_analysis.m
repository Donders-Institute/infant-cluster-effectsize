%% Input of the dataset

clear;

user = getenv('USER');
if isempty(user)
  user = getenv('UserName');
end

switch user
  case 'Didi'    
    data_dir    = 'C:\Users\Didi\Documents\GitHub\Donders Datasets\BeeG dataset\BIDS';
  case 'roboos'
    data_dir    = '/Volumes/Samsung T3/data/di.dcc.DSC_2020.00134_473/bids';
  otherwise
    error('data directory could not be found');
end

addpath(data_dir)
cd(data_dir)

%% Find subject info from BIDS data 

% Read the participants tsv to find subject info
t                   = readtable([data_dir filesep 'participants.tsv'], 'FileType', 'text');
subjectlist         = t.participant_id;

%% Loop over single subjects to do analysis
for ii = 1
%for ii=1:size(subjectlist,1)
    sub = subjectlist{ii};
    do_singlesubject_analysis;
    
    %Display step of analysis
    fprintf('\n')
    disp('------------------')
    disp (['Analysis Subject: ' sub])
    disp('------------------')
    fprintf('\n')
end

%% Now we find and ignore subjects with too many rejected trials

% First we define a trial rejection threshold
threshold         = input('Indicate the threshold for percentage of rejected trials [a number between 0 and 100]');

excluded_participants = [];

for ii=1:size(subjectlist,1)
    sub            = subjectlist{ii};
    % We find the folder containing analysis results for each subject    
    output_dir     = fullfile(fileparts(data_dir), 'results', sub);
    if exist([output_dir filesep 'badtrials.mat'], 'file') && exist([output_dir filesep 'trials.mat'], 'file')
        load([output_dir filesep 'badtrials.mat']); 
        load([output_dir filesep 'trials.mat']);
        rejected_trials = size(badtrials.begsample, 1);
        total_trials    = size(trl_new.begsample, 1);
        percentage_rejected_trials = (rejected_trials/total_trials)*100;
        if percentage_rejected_trials > threshold
            % We have to exclude this participant
            excluded_participants = [excluded_participants, ii];
        end
    else
        % The artefact rejection has not been performed
        warning('Continuing to group analysis but artefact rejection results cannot be found');
    end
end

subjectlist_new = subjectlist;
subjectlist_new(excluded_participants) = [];
 
% we save the exluded participants
results_dir = fullfile(fileparts(data_dir), 'results');
save(fullfile(results_dir, 'excludedparticipants.mat'), 'excluded_participants');


%% Group analysis

do_group_analysis(subjectlist_new, data_dir);

%% Finally we can copy artefact rejection to BIDS, if desired

do_summarize_artefactrejection(subjectlist, data_dir)

