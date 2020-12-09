%% %% Analysis script to perform group analysis of the BeeG dataset

function do_group_analysis(subjectlist)

%% We set the input data directory and the output data directory

user = getenv('USER');
if isempty(user)
  user = getenv('UserName');
end

switch user
  case 'Didi'    
    data_dir    = 'C:\Users\Didi\Documents\GitHub\Donders Datasets\BeeG dataset\results';
    output_dir  = ['C:\Users\Didi\Documents\GitHub\Donders Datasets\BeeG dataset\results' filesep 'group'];    
  case 'roboos'
    data_dir    = '/Volumes/Samsung T3/data/di.dcc.DSC_2020.00134_473/results';
    output_dir  = ['/Volumes/Samsung T3/data/di.dcc.DSC_2020.00134_473/results' filesep 'group'];
  otherwise
    errror('you have to specify the local directories of the data and this code');
end

addpath(data_dir)
cd(data_dir)

if exist(output_dir, 'dir')
    rmdir(output_dir, 's');    
end

mkdir(output_dir);

%% We loop through all subjects and obtain the results of their timelock analysis


% for ii = 1:length(subjectlist)
for ii = 1:5
    folder                  = [data_dir filesep subjectlist{ii}];
    load([folder filesep 'timelock_expected.mat']);
    expected_all(ii)        = { expected }; % We collect all averages in a cell array of structs
    load([folder filesep 'timelock_unexpected.mat']);
    unexpected_all(ii)      = { unexpected };
    load([folder filesep 'timelock_repetition1.mat']);
    repetition1_all(ii)     = { repetition1 };
    load([folder filesep 'timelock_repetition2.mat']);
    repetition2_all(ii)     = { repetition2 };
    load([folder filesep 'timelock_repetition3.mat']);
    repetition3_all(ii)     = { repetition3 };
end

%% calculate grand average for the expected and unexpected stimuli

cfg                      = [];
cfg.channel              = 'all';
cfg.latency              = 'all';
cfg.parameter            = 'avg';
grandavg_expected        = ft_timelockgrandaverage(cfg, expected_all{:});
grandavg_unexpected      = ft_timelockgrandaverage(cfg, unexpected_all{:});

% Then we plot the results
cfg                      = [];
cfg.layout               = 'elec1010.lay';
cfg.interactive          = 'yes';
cfg.showoutline          = 'yes';
cfg.showlabels           = 'yes';
ft_multiplotER(cfg, grandavg_expected, grandavg_unexpected);

% And save the data
save(fullfile(output_dir, 'grandaverage_expected.mat'), 'grandavg_expected');
save(fullfile(output_dir, 'grandaverage_unexpected.mat'), 'grandavg_unexpected');
savefig(gcf, fullfile(output_dir, 'topoplot_grandaverage_expected_unexpected'));


%% Then the same for the repetitions

cfg                      = [];
cfg.channel              = 'all';
cfg.latency              = 'all';
cfg.parameter            = 'avg';
grandavg_repetition1      = ft_timelockgrandaverage(cfg, repetition1_all{:});
grandavg_repetition2      = ft_timelockgrandaverage(cfg, repetition2_all{:});
grandavg_repetition3      = ft_timelockgrandaverage(cfg, repetition3_all{:});

% Then we plot the results
cfg                      = [];
cfg.layout               = 'elec1010.lay';
cfg.interactive          = 'yes';
cfg.showoutline          = 'yes';
cfg.showlabels           = 'yes';
ft_multiplotER(cfg, grandavg_repetition1, grandavg_repetition2, grandavg_repetition3);

% And save the data
save(fullfile(output_dir, 'grandaverage_repetition1.mat'), 'grandavg_repetition1');
save(fullfile(output_dir, 'grandaverage_repetition2.mat'), 'grandavg_repetition2');
save(fullfile(output_dir, 'grandaverage_repetition3.mat'), 'grandavg_repetition3');
savefig(gcf, fullfile(output_dir, 'topoplot_grandaverage_repetitions_expected'));



