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

if ~exist(output_dir, 'dir')
    mkdir(output_dir);  
end



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
cfg.layout                = 'EEG1010.lay';
cfg.highlight             = 'on';
cfg.highlightchannel      = find(stat_expected_unexpected.mask);
cfg.comment               = 'no';
figure; ft_topoplotER(cfg, grandavg_expected)

% And save the data
save(fullfile(output_dir, 'stat_expected_unexpected.mat'), 'stat_expected_unexpected');
savefig(gcf, fullfile(output_dir, 'topoplot_stat_expected_unexpected'));

%% Now we do the statistics with correction for multiple comparisons

% First we need to find neighbouring channels

% This was done with the script below. It is now commented out and we use the struct created from it to continue

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% To start we read in the channels used in the EEG study

% label = expected.cfg.channel;
% 
% % Read sensor positions from the standard 1020 3D format available on fieldtrip
% elec = ft_read_sens('standard_1020.elc', 'senstype', 'eeg');
% 
% if exist('newlayout')
%     clear newlayout
% end
% 
% % Now we loop through all our used channels, to create a new layout struct
% % containing only those channels used in our study
% 
% for ii = 1:length(label)
%       index                     = strcmp(label{ii}, elec.label);
%       newlayout.chanpos(ii, :)  = elec.chanpos(index, :);
%       newlayout.chantype{ii, 1} = elec.chantype{index};
%       newlayout.chanunit{ii, 1} = elec.chanunit{index};
%       newlayout.elecpos(ii, :)  = elec.elecpos(index, :);
%       newlayout.label{ii, 1}    = elec.label{index};
%       newlayout.type            = elec.type;
%       newlayout.unit            = elec.unit;
% end
% 
% % Then we plot this in 3D
% 
% fig = ft_plot_sens(newlayout, 'label', 'label');
% 
% %% Step 2: draw lines in the figure to connect channels and create a neighbours struct
% 
% % We loop through all our channels
% 
% for cc = 1:length(newlayout.label)
%      
%      % find the corresponding label
%      channel_to_connect = newlayout.label{cc};
%          
%      % We then allow the user to provide a cell array of channels to connect
%      % the selected channel to
%      
%      fprintf('\n')
%      disp(['channel to connect = ' channel_to_connect ', channel number ' num2str(cc) ' of ' num2str(length(newlayout.label))]);
%      text = 'provide a row cell array named channels to connect to channel to connect: ';
%      channels = input(text);       
%      
%      % We then create a line between the selected channel and the channels
%      % to connect that channel to
%      
%      connect    = find(strcmp(newlayout.label, channel_to_connect)); 
%      
%      for ci = 1:length(channels)
%          connect1    = find(strcmp(newlayout.label, channels(ci)));
%          linepoints  = [newlayout.chanpos(connect, :); newlayout.chanpos(connect1, :)];
%          line(linepoints(:,1), linepoints(:,2), linepoints(:,3));
%      end
%  
%      % Then we create the neighbours struct
%      
%      selected_neighbours(cc).label       = channel_to_connect;
%      selected_neighbours(cc).neighblabel = channels;
%     
%  end
%  
% % We save the neighbours struct and created image to disk
% save(fullfile(output_dir, 'selected_neighbours.mat'), 'selected_neighbours');
% savefig(gcf, fullfile(output_dir, 'selected_neighbours_plot'));
% 
% % When this step of the code is finished, comment it out to continue with
% % manual optimization in the next step
% 
% %%  Last step: make some manual changes if necessary
% 
% if exist('selected_neighbours.mat')
%     load('selected_neighbours.mat');
% end
% 
% if exist('selected_neighbours_plot.fig')
%     openfig('selected_neighbours_plot.fig');
% end
% 
% % Add a line between two points that were missed:
% 
% connect     = find(strcmp(newlayout.label, 'Fp1'));
% connect1    = find(strcmp(newlayout.label, 'Fp2'));
% linepoints  = [newlayout.chanpos(connect, :); newlayout.chanpos(connect1, :)];
% line(linepoints(:,1), linepoints(:,2), linepoints(:,3));
% 
% % And add the connection to the neighours struct
% selected_neighbours(1).neighblabel  = [selected_neighbours(1).neighblabel, {'Fp2'}];
% selected_neighbours(2).neighblabel  = [selected_neighbours(2).neighblabel, {'Fp1'}];
% 
% % Then we do a sanity check, make sure that neighbours are fully bidirectional
% % I.e. if Cz has FCz as a neighbour, then FCz should have Cz as a neighbour
% 
% for ii = 1:length(selected_neighbours)
%     channel = selected_neighbours(ii).label;
%     for tt  = 1:length(selected_neighbours(ii).neighblabel)
%         idx = find(strcmp(label, selected_neighbours(ii).neighblabel{tt}));
%         if sum(strcmp(selected_neighbours(idx).neighblabel, channel))==0
%           selected_neighbours(idx).neighblabel = [selected_neighbours(idx).neighblabel, {channel}];
%         end
%     end
% end
% 
% % And we replot so that all line are now correct
% 
% close all
% 
% fig = ft_plot_sens(newlayout, 'label', 'label');
% 
% for ii = 1:length(selected_neighbours)
%     channel         = selected_neighbours(ii).label;
%     connect         = find(strcmp(newlayout.label, channel));
%     for tt          = 1:length(selected_neighbours(ii).neighblabel)        
%         connect1    = find(strcmp(newlayout.label, selected_neighbours(ii).neighblabel{tt}));
%         linepoints  = [newlayout.chanpos(connect, :); newlayout.chanpos(connect1, :)];
%         line(linepoints(:,1), linepoints(:,2), linepoints(:,3));
%     end
% end
% 
% % And save this plot and new neighbours
% 
% save(fullfile(output_dir, 'selected_neighbours.mat'), 'selected_neighbours');
% savefig(gcf, fullfile(output_dir, 'selected_neighbours_plot'));

%%  Then we perform the permutation based statistics

if exist([output_dir filesep 'selected_neighbours.mat'], 'file')
    load([output_dir filesep 'selected_neighbours.mat']);
end

cfg                       = [];
cfg.channel               = 'EEG';
cfg.neighbours            = selected_neighbours; % defined as above
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

%% Finally we perform clusterstats for the repetitions of expected stimuli

if exist([output_dir filesep 'selected_neighbours.mat'], 'file')
    load([output_dir filesep 'selected_neighbours.mat']);
end

cfg                       = [];
cfg.channel               = 'EEG';
cfg.neighbours            = selected_neighbours; % defined as above
cfg.parameter             = 'avg';
cfg.method                = 'montecarlo';
cfg.statistic             = 'depsamplesregrT';
cfg.alpha                 = 0.05;
cfg.correctm              = 'cluster';
cfg.correcttail           = 'prob';
cfg.numrandomization      = 1024; % Enlarge this when doing real analysis

%design matrix
Nsub                      = length(subjectlist);
Nrep                      = 3; % three repetitions or "time points" 
design                    = zeros(2, Nrep*Nsub);
design(1, :)              = repmat(1:Nsub, 1, Nrep);  % subject
design(2, :)              = repelem(1:Nrep, Nsub);  % session
cfg.design                = design;
cfg.uvar                  = 1;  % the unit variable: subjects
cfg.ivar                  = 2;  % the independent variable, here number of repetitions

stat_repetitions_clusstats  = ft_timelockstatistics(cfg, repetition1_all{:}, repetition2_all{:}, repetition3_all{:});

% Then plot the result

cfg                       = [];
cfg.style                 = 'blank';
cfg.layout                = 'EEG1010.lay';
cfg.highlight             = 'on';
cfg.highlightchannel      = find(stat_repetitions_clusstats.mask);
cfg.comment               = 'no';
figure; ft_topoplotER(cfg, grandavg_repetition1)

% And save the data
save(fullfile(output_dir, 'stat_repetitions_clusstats.mat'), 'stat_repetitions_clusstats');
savefig(gcf, fullfile(output_dir, 'topoplot_stat_repetitions_clusstats'));


close all
