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


for ii = 1:length(subjectlist)
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
cfg.layout               = 'EEG1010.lay';
cfg.interactive          = 'yes';
cfg.showoutline          = 'yes';
cfg.showlabels           = 'yes';
ft_multiplotER(cfg, grandavg_expected, grandavg_unexpected);

% And save the data
save(fullfile(output_dir, 'grandaverage_expected.mat'), 'grandavg_expected');
save(fullfile(output_dir, 'grandaverage_unexpected.mat'), 'grandavg_unexpected');
savefig(gcf, fullfile(output_dir, 'topoplot_grandaverage_expected_unexpected'));


%% Then the same for the repetitions

cfg                       = [];
cfg.channel               = 'all';
cfg.latency               = 'all';
cfg.parameter             = 'avg';
grandavg_repetition1      = ft_timelockgrandaverage(cfg, repetition1_all{:});
grandavg_repetition2      = ft_timelockgrandaverage(cfg, repetition2_all{:});
grandavg_repetition3      = ft_timelockgrandaverage(cfg, repetition3_all{:});

% Then we plot the results
cfg                       = [];
cfg.layout                = 'EEG1010.lay';
cfg.interactive           = 'yes';
cfg.showoutline           = 'yes';
cfg.showlabels            = 'yes';
ft_multiplotER(cfg, grandavg_repetition1, grandavg_repetition2, grandavg_repetition3);

% And save the data
save(fullfile(output_dir, 'grandaverage_repetition1.mat'), 'grandavg_repetition1');
save(fullfile(output_dir, 'grandaverage_repetition2.mat'), 'grandavg_repetition2');
save(fullfile(output_dir, 'grandaverage_repetition3.mat'), 'grandavg_repetition3');
savefig(gcf, fullfile(output_dir, 'topoplot_grandaverage_repetitions_expected'));

%% Now we perform statistics: permutation based statistics on expected versus unexpected stimuli

% We don't average over time or channels, nor do we choose a specific latency

% First we do not correct for multiple comparisons by clustering

cfg                       = [];
cfg.channel               = 'EEG';
cfg.parameter             = 'avg';
cfg.method                = 'montecarlo';
cfg.statistic             = 'ft_statfun_depsamplesT'; % The samples are dependent (for each subject two conditions)
cfg.alpha                 = 0.05;
cfg.correctm              = 'no'; % No multiple comparisons correction for now
cfg.correcttail           = 'prob'; % distribute the alpha level over both tails by multiplying the probability with a factor two, prior to thresholding it wich cfg.alpha.  
cfg.numrandomization      = 1024; % Enlarge this when doing real analysis

Nsub                      = length(subjectlist);
cfg.design(1,1:2*Nsub)    = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)    = [1:Nsub 1:Nsub];
cfg.ivar                  = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                  = 2; % the 2nd row in cfg.design contains the subject number

stat_expected_unexpected  = ft_timelockstatistics(cfg, expected_all{:}, unexpected_all{:});

% Then plot the result

cfg                       = [];
cfg.style                 = 'blank';
% cfg.layout                = 'EEG1010.lay';
cfg.layout                = 'xyz_layout_1020_33Ch.mat';
cfg.highlight             = 'on';
cfg.highlightchannel      = find(stat_expected_unexpected.mask);
cfg.comment               = 'no';
figure; ft_topoplotER(cfg, grandavg_expected)

% And save the data
save(fullfile(output_dir, 'stat_expected_unexpected.mat'), 'stat_expected_unexpected');
savefig(gcf, fullfile(output_dir, 'topoplot_stat_expected_unexpected'));

%% Now we do the statistics with correction for multiple comparisons

% First we need to find neighbouring channels

cfg                       = [];
cfg.channel               = 'EEG'; 
cfg.method                = 'triangulation';  
cfg.layout                = 'mylayout_32Ch_actiCAP.mat';              
cfg.feedback              = 'yes';                           
neighbours                = ft_prepare_neighbours(cfg, grandavg_expected); % define neighbouring channels

% Then we perform the permutation based statistics as well

cfg                       = [];
cfg.channel               = 'EEG';
cfg.neighbours            = neighbours; % defined as above
cfg.parameter             = 'avg';
cfg.method                = 'montecarlo';
cfg.statistic             = 'ft_statfun_depsamplesT';
cfg.alpha                 = 0.05;
cfg.correctm              = 'cluster';
cfg.correcttail           = 'prob';
cfg.numrandomization      = 1024; % Enlarge this when doing real analysis

Nsub                      = length(subjectlist);
cfg.design(1,1:2*Nsub)    = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)    = [1:Nsub 1:Nsub];
cfg.ivar                  = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                  = 2; % the 2nd row in cfg.design contains the subject number

stat_expected_unexpected_clusstats  = ft_timelockstatistics(cfg, expected_all{:}, unexpected_all{:});

% Then plot the result

cfg                       = [];
cfg.style                 = 'blank';
cfg.layout                = 'EEG1010.lay';
cfg.highlight             = 'on';
cfg.highlightchannel      = find(stat_expected_unexpected_clusstats.mask);
cfg.comment               = 'no';
figure; ft_topoplotER(cfg, grandavg_expected)

% And save the data
save(fullfile(output_dir, 'stat_expected_unexpected_clusstats.mat'), 'stat_expected_unexpected_clusstats');
savefig(gcf, fullfile(output_dir, 'topoplot_stat_expected_unexpected_clusstats'));

close all
