%% %% Analysis script to perform group analysis of the BeeG dataset

do_setpath

% Display step of analysis
fprintf('\n')
disp('------------------------------------')
disp ('Doing group analysis')
disp('------------------------------------')
fprintf('\n')

% Specifying results directory for the group level
output_dir = fullfile(results, 'group');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% First, find and exclude subjects if too many trials had to be rejected

% First we define a trial rejection threshold
threshold = input('Indicate the threshold for percentage of rejected trials [a number between 0 and 100]');

excluded_participants = [];

for ii = 1:size(subjectlist,1)
    % We find the folder containing analysis results for each subject, those are the input for the group analysis
    sub            = subjectlist{ii};
    input_dir     = fullfile(results, sub);
    if exist([input_dir filesep 'badtrials.mat'], 'file') && exist([input_dir filesep 'trials.mat'], 'file')
        load([input_dir filesep 'badtrials.mat']);
        load([input_dir filesep 'trials.mat']);
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
save(fullfile(output_dir, 'excludedparticipants.mat'), 'excluded_participants');

%% We loop through all subjects and obtain the results of their timelock analysis
load(fullfile(output_dir, 'excludedparticipants.mat'), 'excluded_participants');

for ii = 1:length(subjectlist_new)
    folder                  = [results filesep subjectlist_new{ii}];
    load([folder filesep 'timelock_expected.mat']);
    expected_all(ii)        = { expected }; % We collect all averages in a cell array of structs
    load([folder filesep 'timelock_unexpected.mat']);
    unexpected_all(ii)      = { unexpected };
%     load([folder filesep 'timelock_repetition1.mat']);
%     repetition1_all(ii)     = { repetition1 };
%     load([folder filesep 'timelock_repetition2.mat']);
%     repetition2_all(ii)     = { repetition2 };
%     load([folder filesep 'timelock_repetition3.mat']);
%     repetition3_all(ii)     = { repetition3 };
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

Nsub                      = length(subjectlist_new);
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

%% Perform cluster-based permutation statistics

load(fullfile(scripts, 'selected_neighbours.mat'));

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

Nsub                      = length(subjectlist_new);
cfg.design(1,1:2*Nsub)    = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)    = [1:Nsub 1:Nsub];
cfg.ivar                  = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                  = 2; % the 2nd row in cfg.design contains the subject number

stat_expected_unexpected_clusstats  = ft_timelockstatistics(cfg, expected_all{:}, unexpected_all{:});

%% Plot the results
% For this purpose, first calculate the difference between conditions
cfg                       = [];
cfg.operation             = 'subtract';
cfg.parameter             = 'avg';
grandavg_diff_exp_vs_unexp= ft_math(cfg, grandavg_expected, grandavg_unexpected);

% Find clusters with a 5% two-sided cutoff based on the cluster p-values
pos_cluster_pvals = [stat_expected_unexpected_clusstats.posclusters(:).prob];
pos_clust = find(pos_cluster_pvals < 0.025);
pos = ismember(stat_expected_unexpected_clusstats.posclusterslabelmat, pos_clust);

neg_cluster_pvals = [stat_expected_unexpected_clusstats.negclusters(:).prob];
neg_clust = find(neg_cluster_pvals < 0.025);
neg = ismember(stat_expected_unexpected_clusstats.negclusterslabelmat, neg_clust);

% % Alternatively, plot only the first positive/negative cluster
% pos = stat_expected_unexpected_clusstats.posclusterslabelmat ==1;
% neg = stat_expected_unexpected_clusstats.negclusterslabelmat ==1;

% Set plotting specifications
timestep = 0.05; % plot every 0.05 sec intervals
sampling_rate = 500; % set sampling frequency %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTE: make dynamic
sample_count = length(stat_expected_unexpected_clusstats.time);
j = [stat_expected_unexpected_clusstats.time(1):timestep:stat_expected_unexpected_clusstats.time(end)];% start of each interval for plotting in seconds
m = [1:timestep*sampling_rate:sample_count]; % start of each interval for plotting in sample points

[i1,i2] = match_str(grandavg_diff_exp_vs_unexp.label, stat_expected_unexpected_clusstats.label); % matching channel labels

for k = 1:30
    subplot(6,5,k);
    cfg = [];
    cfg.xlim = [j(k) j(k+1)]; % current interval
    cfg.zlim = [-10 10]; % set minimum and maximum z-axis
    pos_int = zeros(numel(grandavg_diff_exp_vs_unexp.label),1);
    neg_int = zeros(numel(grandavg_diff_exp_vs_unexp.label),1);
    pos_int(i1) = all(pos(i2, m(k):m(k+1)),2); % determine which channel are in a cluster throughout the current time interval (pos cluster)
    neg_int(i1) = all(neg(i2, m(k):m(k+1)),2); % determine which channel are in a cluster throughout the current time interval (neg cluster)
    
    cfg.highlight = 'on';
    cfg.highlightchannel = find(pos_int | neg_int); % highlight channels belonging to a cluster
    cfg.highlightcolor     = [1 1 1]; % highlight marker color (default = [0 0 0] (black))
    cfg.comment = 'xlim';
    cfg.commentpos = 'title';
    cfg.layout =  'EEG1010.lay';
    cfg.interactive = 'no';
    ft_topoplotER(cfg, grandavg_diff_exp_vs_unexp)
    
end

% And save the data
save(fullfile(output_dir, 'stat_expected_unexpected_clusstats.mat'), 'stat_expected_unexpected_clusstats');
savefig(gcf, fullfile(output_dir, 'topoplot_stat_expected_unexpected_clusstats'));

%% Determine effect size
% Calculate Cohen's d for the average difference in the respective cluster and determine at which channel/time this Cohen's d
% is maximal 

% First for the positive cluster
effect_window_pos = stat_expected_unexpected_clusstats.posclusterslabelmat==1;

for ii = 1:size(subjectlist_new,1)
    folder                  = [results filesep subjectlist_new{ii}];
    load([folder filesep 'timelock_expected.mat']);
    load([folder filesep 'timelock_unexpected.mat']);
    
    a = expected.avg(effect_window_pos);
    b = unexpected.avg(effect_window_pos);
    c = a-b;
    Pos.ERP_Diff_alltimechan(ii,:) =c;
    Pos.ERP_Diff(ii) = nanmean([a - b]);
    clear expected unexpected a b c

end
Pos.stdev_ERP_diff = std(Pos.ERP_Diff);
Pos.mean_ERP_diff = mean(Pos.ERP_Diff);
Pos.cohensd_ERP_diff = Pos.mean_ERP_diff/Pos.stdev_ERP_diff;

for t = 1:size(Pos.ERP_Diff_alltimechan,2)   
    Pos.cohensd_ERP_Diff_alltimechan(t)=(nanmean(Pos.ERP_Diff_alltimechan(:,t)))/(std(Pos.ERP_Diff_alltimechan(:,t)));
end
[ Pos.cohensd_ERP_Diff_Max, idx]= max(Pos.cohensd_ERP_Diff_alltimechan);

[Pos.row,Pos.col] = find(stat_expected_unexpected_clusstats.posclusterslabelmat==1);
 
folder                  = [results filesep subjectlist_new{1}]; % load an example participant's ERP for determining channel and time information
load([folder filesep 'timelock_expected.mat']);

fprintf('\n')
disp('~~~~~')
disp(['Contrast: Standard vs. Oddball: ']);
disp(['Effect size Cohens d of average positive cluster is ' num2str(Pos.cohensd_ERP_diff)])
disp(['Maximum effect size is ' num2str(Pos.cohensd_ERP_Diff_Max) ' at channel ' expected.label{Pos.row(idx)} ' and at time ' num2str(expected.time(Pos.col(idx))) ' sec']);
disp('~~~~~')
fprintf('\n')

save(fullfile(output_dir, 'EffectSize_Pos.mat'), 'Pos');

% Then for the negative clustesr
effect_window_neg = stat_expected_unexpected_clusstats.negclusterslabelmat==1;


for ii = 1:size(subjectlist_new,1)
    folder                  = [results filesep subjectlist_new{ii}];
    load([folder filesep 'timelock_expected.mat']);
    load([folder filesep 'timelock_unexpected.mat']);
    
    a = expected.avg(effect_window_neg);
    b = unexpected.avg(effect_window_neg);
    c = a-b;
    Neg.ERP_Diff_alltimechan(ii,:) =c;
    Neg.ERP_Diff(ii) = nanmean([a - b]);
    clear expected unexpected a b c

end
Neg.stdev_ERP_diff = std(Neg.ERP_Diff);
Neg.mean_ERP_diff = mean(Neg.ERP_Diff);
Neg.cohensd_ERP_diff = Neg.mean_ERP_diff/Neg.stdev_ERP_diff;

for t = 1:size(Neg.ERP_Diff_alltimechan,2)   
    Neg.cohensd_ERP_Diff_alltimechan(t)=(nanmean(Neg.ERP_Diff_alltimechan(:,t)))/(std(Neg.ERP_Diff_alltimechan(:,t)));
end
[ Neg.cohensd_ERP_Diff_Max, idx]= min(Neg.cohensd_ERP_Diff_alltimechan);

[Neg.row,Neg.col] = find(stat_expected_unexpected_clusstats.negclusterslabelmat==1);
 
folder                  = [results filesep subjectlist_new{1}]; % load an example participant's ERP for determining channel and time information
load([folder filesep 'timelock_expected.mat']);

fprintf('\n')
disp('~~~~~')
disp(['Contrast: Standard vs. Oddball: '])
disp(['Effect size Cohens d of average positive cluster is ' num2str(Neg.cohensd_ERP_diff)])
disp(['Maximum effect size is ' num2str(Neg.cohensd_ERP_Diff_Max) ' at channel ' expected.label{Neg.row(idx)} ' and at time ' num2str(expected.time(Neg.col(idx))) ' sec']);
disp('~~~~~')
fprintf('\n')

save(fullfile(output_dir, 'EffectSize_Neg.mat'), 'Neg');

%% Plot ERP timecourse of the electrode with the maximal effect size

% Determine variability between participants
cfg                      = [];
cfg.channel              = 'all';
cfg.latency              = 'all';
cfg.parameter            = 'avg';
cfg.keepindividual       = 'yes';
grandavg_expected_all        = ft_timelockgrandaverage(cfg, expected_all{:});
grandavg_unexpected_all     = ft_timelockgrandaverage(cfg, unexpected_all{:});

se_grandavg_expected = squeeze(nanstd(grandavg_expected_all.individual/sqrt(length(subjectlist_new))));
se_grandavg_unexpected = squeeze(nanstd(grandavg_unexpected_all.individual/sqrt(length(subjectlist_new))));

% Set color specifications
colour_code = {'b','g', 'm', 'r', 'c', 'k', ':k', 'g', 'm'};
shaded_area = {[0, 0, 1],[0, 1, 0],[1 0 1], [1 0 0], [0 1 1], [0, 1, 1],[0, 1, 0],[0, 0, 0],[0, 0, 0]};

% Plot 
figure;
% ERP of channel with maximum effect size of positive Cluster
subplot(1,2,1)
% Condition 1
plot(grandavg_expected.time,grandavg_expected.avg(Pos.row(idx),:),colour_code{1}, 'LineWidth', 1.5)

mean_expected = grandavg_expected.avg(Pos.row(idx),:);
se_expected = se_grandavg_expected(Pos.row(idx),:);
patch([grandavg_expected.time, fliplr(grandavg_expected.time)], [mean_expected-se_expected, fliplr(mean_expected+se_expected)],  shaded_area{1}, 'edgecolor', 'none', 'FaceAlpha', .3);

hold all
% Condition 2
plot(grandavg_unexpected.time,grandavg_unexpected.avg(Pos.row(idx),:),colour_code{4}, 'LineWidth', 1.5)

mean_unexpected = grandavg_unexpected.avg(Pos.row(idx),:);
se_unexpected = se_grandavg_unexpected(Pos.row(idx),:);
patch([grandavg_expected.time, fliplr(grandavg_expected.time)], [mean_unexpected-se_unexpected, fliplr(mean_unexpected+se_unexpected)],  shaded_area{4}, 'edgecolor', 'none', 'FaceAlpha', .3);

xlabel('Time [s]');
ylabel('Amplitude [mV]');
ylim([-15 15])
line([-.5 1],[ 0 0], 'Color', [0 0 0],'LineStyle', ':')
title(['Maximum effect of positive cluster at channel ' grandavg_expected.label(Pos.row(idx))])

% ERP of channel with maximum effect size of positive Cluster
subplot(1,2,2)
% Condition 1
plot(grandavg_expected.time,grandavg_expected.avg(Neg.row(idx),:),colour_code{1}, 'LineWidth', 1.5)

mean_expected = grandavg_expected.avg(Neg.row(idx),:);
se_expected = se_grandavg_expected(Neg.row(idx),:);
patch([grandavg_expected.time, fliplr(grandavg_expected.time)], [mean_expected-se_expected, fliplr(mean_expected+se_expected)],  shaded_area{1}, 'edgecolor', 'none', 'FaceAlpha', .3);

hold all
% Condition 2
plot(grandavg_unexpected.time,grandavg_unexpected.avg(Neg.row(idx),:),colour_code{4}, 'LineWidth', 1.5)

mean_unexpected = grandavg_unexpected.avg(Neg.row(idx),:);
se_unexpected = se_grandavg_unexpected(Neg.row(idx),:);
patch([grandavg_expected.time, fliplr(grandavg_expected.time)], [mean_unexpected-se_unexpected, fliplr(mean_unexpected+se_unexpected)],  shaded_area{4}, 'edgecolor', 'none', 'FaceAlpha', .3);

xlabel('Time [s]');
ylabel('Amplitude [mV]');
ylim([-15 15])
line([-.5 1],[ 0 0], 'Color', [0 0 0],'LineStyle', ':')
title(['Maximum effect of negative cluster at channel ' grandavg_expected.label(Neg.row(idx))])

%% Finally we perform clusterstats for the repetitions of expected stimuli

if exist([output_dir filesep 'selected_neighbours.mat'], 'file')
    load([output_dir filesep 'selected_neighbours.mat']);
end

cfg                       = [];
cfg.channel               = 'EEG';
cfg.neighbours            = selected_neighbours; % defined as above
cfg.parameter             = 'avg';
cfg.method                = 'montecarlo';
cfg.statistic             = 'depsamplesregrT'; % Perhaps Spearman would be better here
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
