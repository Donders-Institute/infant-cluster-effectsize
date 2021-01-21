%% Copy rejected trials and rejected channels to BIDS

% This script should be run when the artefact rejection is complete. It will
% incorporate the rejected channels and trials in a derivatives subfolder.
%
% Specifically, number of rejected trials and excluded participants will go in a
% participants tsv
%
% Begsample, endsample, and round of artefact rejection of each rejected trial go in
% a badtrials tsv file
%
% Rejected channels go in a channels tsv

% We create a derivatives directory
if ~exist(derivatives, 'dir')
    mkdir(derivatives);
end

% We find the excluded participants
results_dir = fullfile(fileparts(bidsroot), 'results');
if exist([results_dir filesep 'excludedparticipants.mat'], 'file')
    load([results_dir filesep 'excludedparticipants.mat']);
else
    error('number of exluded participants could not be found');
end

% We initiate the exluded trials vector (to be filled out later)
excluded_trials                   = nan(length(subjectlist), 1);

% We loop over subjects
for ii=1:size(subjectlist,1)
    sub = subjectlist{ii};
    
    % and find the output dir that contains the stored badchannels and
    % badtrials
    
    output_dir = fullfile(fileparts(bidsroot), 'results', sub);
    % From here we load the badtrials and badchannels
    if exist([output_dir filesep 'badtrials.mat'], 'file')
        load([output_dir filesep 'badtrials.mat']);
        % The we read the participants tsv, to which we should add
        % the rejected trial info
        participants = read_tsv([bidsroot filesep 'participants.tsv']);
        
        % We write the included column
        included = repmat({'Yes'}, length(subjectlist), 1);
        included(excluded_participants, :)  = {'No'};
        participants.included = included;
        
        excluded_trials(ii, 1)            = size(badtrials.begsample, 1);
        participants.excluded_trials      = excluded_trials;
        
        % Then we update the current tsv file
        write_tsv([derivatives filesep 'participants.tsv'], participants);
        
        % We also have to update the corresponding json file
        participants_json      = read_json([bidsroot filesep 'participants.json']);
        if ~isfield(participants_json, 'included')
            % the field does not yet exist, we create it
            participants_json.included.description = 'Whether participant has been included in the final analysis based on percentage of rejected trials';
            participants_json.included.levels      = {'Yes: participant was included in the analysis',...
                'No: participant was excluded from the analysis'};
            % and re-write the json file
            write_json([derivatives filesep 'participants.json'], participants_json);
        end
        if ~isfield(participants_json, 'excluded_trials')
            % the field does not yet exist, we create it
            participants_json.excluded_trials.description = 'number of rejected trials during analysis';
            % and re-write the json file
            write_json([derivatives filesep 'participants.json'], participants_json);
        end
        
        % Secondly, we create an events tsv file containing begsample,
        % endsample, and rejection round of each rejected trial
        
        % Only if badtrials is not empty
        
        if ~isempty(badtrials)
            
            if exist([derivatives filesep sub filesep 'eeg' filesep sub '_badtrials.tsv'], 'file')
                % The file already exists: we update it
                write_tsv([derivatives filesep sub filesep 'eeg' filesep sub '_badtrials.tsv'], badtrials);
            else
                % We have to make the folder and file
                mkdir([derivatives filesep sub filesep 'eeg'])
                write_tsv([derivatives filesep sub filesep 'eeg' filesep sub '_badtrials.tsv'], badtrials);
                
                % Then we also make a corresponding json file
                badtrials_json.begsample.description       = 'Sample where event begins (measured from start of recording)';
                badtrials_json.begsample.units             = 'sample number';
                badtrials_json.endsample.description       = 'Sample where event ends (measured from start of recording)';
                badtrials_json.endsample.units             = 'sample number';
                badtrials_json.rejection_round.description = 'Round of the artefact rejection process where the artefact has been detected';
                badtrials_json.rejection_round.levels      = {'Visual rejection 1: the first, rough visual rejection round',...
                    'Visual rejection 2: the last, fine-tuning rejection round'};
                
                % Then we write json file
                write_json([derivatives filesep sub filesep 'eeg' filesep sub '_badtrials.json'], badtrials_json);
            end
        end
        
    else
        % Apparently no artefact rejection has been performed yet
        str = ['no artefact rejection has been performed yet for subject ' sub];
        warning(str)
    end
    
    % Now for the rejected channels
    if exist([output_dir filesep 'badchannelnumbers.mat'], 'file')
        load([output_dir filesep 'badchannelnumbers.mat']);
        
        % We check if the folder has already been created
        if ~exist([derivatives filesep sub filesep 'eeg'])
            mkdir([derivatives filesep sub filesep 'eeg']);
        end
        
        channels_folder         = dir([bidsroot filesep sub filesep 'eeg' filesep '*channels*.tsv']);
        channels_name           = [bidsroot filesep sub filesep 'eeg' filesep channels_folder.name];
        channels_tsv            = read_tsv(channels_name);
        
        % We add a status column to the channels tsv or replace the existing one
        % with rejected channels
        
        status                  = repmat({'Good'}, length(channels_tsv.name), 1);
        status(badchannels, :)  = {'Bad'};
        channels_tsv.status     = status;
        
        % Then we update the current tsv file
        write_tsv([derivatives filesep sub filesep 'eeg' filesep sub '_channels.tsv'], channels_tsv);
    else
        % do nothing
    end
end

