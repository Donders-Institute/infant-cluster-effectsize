%% Do the complete analysis
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
% Although this script can be executed as it is, it is more likely that the analysis
% will be done over multiple days and that the do_singlesubject_analysis script will
% be called repeatedly. This script mainly serves to show how the different scripts
% fit together.

clear

do_setpath

%% Read the subject info from the BIDS dataset
t           = readtable([bidsroot filesep 'participants.tsv'], 'FileType', 'text');
subjectlist = t.participant_id;

%% Loop over single subjects to do the analysis
for ii = 1:size(subjectlist,1)
  sub = subjectlist{ii};
  do_singlesubject_analysis;
end

%% Do the group analysis

do_group_analysis;
