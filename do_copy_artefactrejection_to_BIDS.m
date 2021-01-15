%% Copy rejected trials and rejected channels to BIDS

function do_copy_artefactrejection_to_BIDS(subjectlist, data_dir)
% This script should be run when the artefact rejection is complete. 
% It will incorporate the rejected channels and trials in the BIDS
% structure

% Specifically, number of rejected trials will go in the participants tsv
% Begsample, endsample, and round of artefact rejection of each rejected
% trial go in the derivatives folder in a badtrials tsv file

% Rejected channels go in the channels tsv

% We loop over subjects
for ii=1:size(subjectlist,1)
    sub = subjectlist{ii};
    % and find the output dir that contains the stored badchannels and
    % badtrials
    output_dir = fullfile(fileparts(data_dir), 'results', sub);
    % From here we load the badtrials and badchannels
    if exist([output_dir filesep 'badtrials.mat'], 'file')
        load([output_dir filesep 'badtrials.mat']); 
        % The we read the participants tsv, to which we should add
        % the rejected trial info
        participants = read_tsv([data_dir filesep 'participants.tsv']);
    
        % Check if the trials column already exists, and add excluded trial number
        if sum(ismember(participants.Properties.VariableNames, 'excluded_trials')) == 0
            % The column does not yet exist: we create it
            excluded_trials                   = nan(length(subjectlist), 1);
            excluded_trials(ii, 1)            = length(badtrials.begsample); 
            participants.excluded_trials      = excluded_trials;
        else
            % The column exists: we add our relevant value to it
            participants.excluded_trials(ii,1) = length(badtrials.begsample);
        end

        % Then we update the current tsv file
        write_tsv([data_dir filesep 'participants.tsv'], participants);
    
        % We also have to update the corresponding json file
        participants_json      = read_json([data_dir filesep 'participants.json']);
        if ~isfield(participants_json, 'excluded_trials')
            % the field does not yet exist, we create it
            participants_json.excluded_trials.description = 'number of rejected trials during analysis';
            % and re-write the json file
            write_json([data_dir filesep 'participants.json'], participants_json);
        end
        
        % Secondly, we create an events tsv file containing begsample,
        % endsample, and rejection round of each rejected trial
        
        if exist([data_dir filesep 'derivatives' filesep sub filesep 'eeg' filesep sub '_badtrials.tsv'], 'file')
            % The file already exists: we update it
            write_tsv([data_dir filesep 'derivatives' filesep sub filesep 'eeg' filesep sub '_badtrials.tsv'], badtrials);
        else
            % We have to make the folder and file
            mkdir([data_dir filesep 'derivatives' filesep sub filesep 'eeg'])
            write_tsv([data_dir filesep 'derivatives' filesep sub filesep 'eeg' filesep sub '_badtrials.tsv'], badtrials);
            
            % Then we also make a corresponding json file
            badtrials_json.begsample.description       = 'Sample where event begins (measured from start of recording)';
            badtrials_json.begsample.units             = 'sample number';
            badtrials_json.endsample.description       = 'Sample where event ends (measured from start of recording)';
            badtrials_json.endsample.units             = 'sample number';
            badtrials_json.rejection_round.description = 'Round of the artefact rejection process where the artefact has been detected';
            badtrials_json.rejection_round.levels      = {'Visual rejection 1: the first, rough visual rejection round',...
                                                          'Visual rejection 2: the last, fine-tuning rejection round'};
                                                      
            % Then we write json file
            write_json([data_dir filesep 'derivatives' filesep sub filesep 'eeg' filesep sub '_badtrials.json'], badtrials_json);
        end    
        
    else
        % Apparently no artefact rejection has been performed yet
        str = ['no artefact rejection has been performed yet for subject ' sub];
        warning(str)    
    end
    
    % Now for the rejected channels
    if exist([output_dir filesep 'badchannelnumbers.mat'], 'file')
        load([output_dir filesep 'badchannelnumbers.mat']);
        
        channels_folder         = dir([data_dir filesep sub filesep 'eeg' filesep '*channels*.tsv']);
        channels_name           = [data_dir filesep sub filesep 'eeg' filesep channels_folder.name];
        channels_tsv            = read_tsv(channels_name);

        % We add a status column to the channels tsv or replace the existing one
        % with rejected channels

        status                  = repmat({'Good'}, length(channels_tsv.name), 1);
        status(badchannels, :)  = {'Bad'};  
        channels_tsv.status     = status;

        % Then we update the current tsv file
        write_tsv(channels_name, channels_tsv);
    else
        % do nothing            
    end       
end


    

    