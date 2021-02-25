%% Input of the dataset

clear

do_setpath

%% Find subject info from BIDS data

cd(bidsroot)

% Read the participants tsv to find subject info
t           = readtable([bidsroot filesep 'participants.tsv'], 'FileType', 'text');
subjectlist = t.participant_id;

% Make a selection to speed testing up
subjectlist = subjectlist(41:end);

%% Loop over single subjects to do analysis
for ii = 1:size(subjectlist,1)
    sub = subjectlist{ii};
    do_singlesubject_analysis;
end

%% Group analysis

do_group_analysis;

%% Finally we can copy artefact rejection to BIDS, if desired

do_summarize_artefactrejection;
