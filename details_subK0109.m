% this script contains some subject specific variables and details.
% Most of them are (almost) the same over subjects, but if subjects were
% inconsistent, this is the place to deal with it.

subject = 'K0109';
subjectnumber = 1;
sub     = ['sub-' subject];

% this is where the intermediate results will be written
% data_dir points to the BIDS dataset
output_dir = fullfile(fileparts(data_dir), 'results', sub);
data_dir   = data_dir;

% The rejected trials and rejected channels can be filled out here:
badchannels = {};
badtrials   = [];

% If nothing has been defined by the user (badchannels and badtrials are empty)
% we read info from the BIDS data structure or output_dir

if isempty(badtrials)
    % none have been defined, we try to read them from the BIDS and output
    % data structure
    
    % The trial info should be in the participants tsv
    participants = read_tsv([data_dir filesep 'participants.tsv']);

    if ~isnan(participants.excluded_trials(subjectnumber, 1))  
        % The number of rejected trials was stored in the participants tsv
        number_trials_rejected = participants.excluded_trials(subjectnumber, 1);
        
        % A table containing rejected trials begsample, endsample, and
        % round of rejection is stored in the output directory
        
        if exist([output_dir filesep 'badtrials.mat'], 'file')
            badtrials = load([output_dir filesep 'badtrials.mat']); 
        end    
    end
else
    % we leave the data as filled out by the user    
end

if isempty(badchannels)
    % none have been defined, we try to read them from the BIDS and output
    % data structure
    % Now for the bad channels, they are in the channels tsv
    channels_folder   = dir([data_dir filesep sub filesep 'eeg' filesep '*channels*.tsv']);
    channels_name     = [data_dir filesep sub filesep 'eeg' filesep channels_folder.name];
    channels_tsv      = read_tsv(channels_name);

    if sum(ismember(channels_tsv.Properties.VariableNames, 'status'))~=0
        % the channels have been analyzed
        if exist([output_dir filesep 'badchannels.mat'], 'file')
            badchannels = load([output_dir filesep 'badchannels.mat']); 
        end    
    end
else
    % we leave the data as filled out by the user  
end

% Finally we ask whether or not to redo artefact rejection for this subject
yn = input('Do you want to redo artefact rejection? [press y / n]','s');
if strcmp(yn,'y')==1
    do_artefact_rejection = 1;
elseif strcmp(yn,'n')==1
    do_artefact_rejection = 0;
end
