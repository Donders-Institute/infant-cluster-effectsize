

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
    errror('you have to specify the local directories of the data and this code');
end

addpath(data_dir)
cd(data_dir)

%% Find subject info from BIDS data 

% Read the participants tsv to find subject info
t                   = readtable([data_dir filesep 'participants.tsv'], 'FileType', 'text');
subjectlist         = t.participant_id;

%% Loop over single subjects to do analysis

% for ii=1:size(subjectlist,1)
for ii = 1
    sub = subjectlist{ii};
    do_singlesubject_analysis(sub);
    
    %Display step of analysis
    fprintf('\n')
    disp('------------------')
    disp (['Analysis Subject: ' sub])
    disp('------------------')
    fprintf('\n')
end

%% Group analysis

do_group_analysis(subjectlist);
