%% %% Analysis script to perform group analysis of the example infant EEG dataset

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

%% First, find and exclude subjects for who too many trials had to be rejected

% Define a trial rejection threshold
threshold = input('Indicate the threshold for percentage of rejected trials [a number between 0 and 100]');

excluded_participants = [];

for ii = 1:size(subjectlist,1)
    sub            = subjectlist{ii};
    input_dir     = fullfile(results, sub);
    
    % Find information on how many trials were rejected & calculate percent
    % of rejected trials
    if exist([input_dir filesep 'badtrials.mat'], 'file') && exist([input_dir filesep 'trials.mat'], 'file')
        load([input_dir filesep 'badtrials.mat']);
        load([input_dir filesep 'trials.mat']);
        rejected_trials = size(badtrials.begsample, 1);
        total_trials    = size(trl_new.begsample, 1);
        percentage_rejected_trials = (rejected_trials/total_trials)*100;
        % Exclude participants if more than the defined threshold of trials
        % were rejected during artifact rejection
        if percentage_rejected_trials > threshold
            excluded_participants = [excluded_participants, ii];
        end
    else
        % Give warning if artefact rejection has not been performed yet
        warning('Continuing to group analysis but artefact rejection results cannot be found');
    end
end

% create updated subject list including only those participants with
% sufficient artifact-free data
subjectlist_new = subjectlist;
subjectlist_new(excluded_participants) = [];

% Save which participants were excluded
save(fullfile(output_dir, 'excludedparticipants.mat'), 'excluded_participants');
save(fullfile(output_dir, 'subjectlist_new.mat'), 'subjectlist_new');

%% Calculate grand average ERP
% Do so by averaging time-locked data across participants
load(fullfile(output_dir, 'subjectlist_new.mat'), 'subjectlist_new');

for ii = 1:length(subjectlist_new)
    folder                  = [results filesep subjectlist_new{ii}];
    load([folder filesep 'timelock_standard.mat']);
    standard_all(ii)        = { standard };
    load([folder filesep 'timelock_oddball.mat']);
    oddball_all(ii)      = { oddball };
end

% Calculate grand average for both conditions (standard, oddball)
% separately

cfg                      = [];
cfg.channel              = 'all';
cfg.latency              = 'all';
cfg.parameter            = 'avg';
grandavg_standard        = ft_timelockgrandaverage(cfg, standard_all{:});
grandavg_oddball         = ft_timelockgrandaverage(cfg, oddball_all{:});

% Plot the results
cfg                      = [];
cfg.layout               = 'EEG1010.lay';
cfg.interactive          = 'yes';
cfg.showoutline          = 'yes';
cfg.showlabels           = 'yes';
ft_multiplotER(cfg, grandavg_standard, grandavg_oddball);

% Save the data
save(fullfile(output_dir, 'grandaverage_standard.mat'), 'grandavg_standard');
save(fullfile(output_dir, 'grandaverage_oddball.mat'), 'grandavg_oddball');
savefig(gcf, fullfile(output_dir, 'topoplot_grandaverage_standard_oddball'));

%% Perform cluster-based permutation statistics correcting for multiple comparisons

load(fullfile(scripts, 'selected_neighbours.mat'));

cfg                       = [];
cfg.channel               = 'EEG';
cfg.neighbours            = selected_neighbours; % as defined for this channel layout
cfg.parameter             = 'avg';
cfg.method                = 'montecarlo';
cfg.statistic             = 'ft_statfun_depsamplesT';
cfg.alpha                 = 0.05;
cfg.correctm              = 'cluster';
cfg.correcttail           = 'prob';
cfg.numrandomization      = 1024; 

Nsub                      = length(subjectlist_new);
cfg.design(1,1:2*Nsub)    = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)    = [1:Nsub 1:Nsub];
cfg.ivar                  = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                  = 2; % the 2nd row in cfg.design contains the subject number

stat_standard_oddball_clusstats  = ft_timelockstatistics(cfg, standard_all{:}, oddball_all{:});

% Save the data
save(fullfile(output_dir, 'stat_standard_oddball_clusstats.mat'), 'stat_standard_oddball_clusstats');

%% Plot the results of the cluster-based permutation test
load(fullfile(output_dir, 'stat_standard_oddball_clusstats.mat'), 'stat_standard_oddball_clusstats');

% Plot displaying t- and p-value distribution across channels and time
plot_clus = zeros(size(stat_standard_oddball_clusstats.prob));
plot_clus(stat_standard_oddball_clusstats.negclusterslabelmat ==1) = -1;
plot_clus(stat_standard_oddball_clusstats.posclusterslabelmat ==1) = 1;

figure
subplot(2,1,1)
imagesc(stat_standard_oddball_clusstats.time, 1:size(stat_standard_oddball_clusstats.label,1),plot_clus)
colormap(jet)
colorbar

title(['Largest positive and negative cluster']);
subplot(2,1,2)
imagesc(stat_standard_oddball_clusstats.time, 1:size(stat_standard_oddball_clusstats.label,1),  stat_standard_oddball_clusstats.stat)
colorbar
title(['T-values per channel x time']);
savefig(gcf, fullfile(output_dir, 'T_and_Pvalues_stat_standard_oddball_clusstats'));

% Plot displaying topographic maps across time bins highlighting channel/time as part of clusters

% For this purpose, calculate the difference between conditions
cfg                       = [];
cfg.operation             = 'subtract';
cfg.parameter             = 'avg';
grandavg_diff_standard_oddball= ft_math(cfg, grandavg_standard, grandavg_oddball);

% Find clusters with a 5% two-sided cutoff based on the cluster p-values
pos_cluster_pvals       = [stat_standard_oddball_clusstats.posclusters(:).prob];
pos_clust               = find(pos_cluster_pvals < 0.025);
pos                     = ismember(stat_standard_oddball_clusstats.posclusterslabelmat, pos_clust);

neg_cluster_pvals       = [stat_standard_oddball_clusstats.negclusters(:).prob];
neg_clust               = find(neg_cluster_pvals < 0.025);
neg                     = ismember(stat_standard_oddball_clusstats.negclusterslabelmat, neg_clust);

% % Alternatively, plot only the first positive/negative cluster
% pos = stat_expected_unexpected_clusstats.posclusterslabelmat ==1;
% neg = stat_expected_unexpected_clusstats.negclusterslabelmat ==1;

% Set plotting specifications
timestep                = 0.05; % plot every 0.05 sec intervals
sampling_rate           = 500; % set sampling frequency 
sample_count            = length(stat_standard_oddball_clusstats.time);
j                       = [stat_standard_oddball_clusstats.time(1):timestep:stat_standard_oddball_clusstats.time(end)]; % start of each interval for plotting in seconds
m                       = [1:timestep*sampling_rate:sample_count]; % start of each interval for plotting in sample points

[i1,i2] = match_str(grandavg_diff_standard_oddball.label, stat_standard_oddball_clusstats.label); % matching channel labels

figure
for k = 1:30
    subplot(6,5,k);
    cfg                 = [];
    cfg.xlim            = [j(k) j(k+1)]; % current interval
    cfg.zlim            = [-10 10]; % set minimum and maximum z-axis
    pos_int             = zeros(numel(grandavg_diff_standard_oddball.label),1);
    neg_int             = zeros(numel(grandavg_diff_standard_oddball.label),1);
    pos_int(i1)         = all(pos(i2, m(k):m(k+1)),2); % determine which channels are in a cluster throughout the current time interval (pos cluster)
    neg_int(i1)         = all(neg(i2, m(k):m(k+1)),2); % determine which channels are in a cluster throughout the current time interval (neg cluster)
    
    cfg.highlight       = 'on';
    cfg.highlightchannel= find(pos_int | neg_int); % highlight channels belonging to a cluster
    cfg.highlightcolor  = [1 1 1]; % highlight marker color (default = [0 0 0] (black))
    cfg.comment         = 'xlim';
    cfg.commentpos      = 'title';
    cfg.layout          =  'EEG1010.lay';
    cfg.interactive     = 'no';
    ft_topoplotER(cfg, grandavg_diff_standard_oddball)
end

% Save the figure
savefig(gcf, fullfile(output_dir, 'topoplot_stat_standard_oddball_clusstats'));

%% Determine effect size

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Option 1: Calculate Cohen's d for the average difference in the respective cluster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First for the positive cluster
effect_window_pos = stat_standard_oddball_clusstats.posclusterslabelmat==1;

% Calculate pairwise difference between conditions for each participant
for ii = 1:size(subjectlist_new,1)
    folder                  = [results filesep subjectlist_new{ii}];
    load([folder filesep 'timelock_standard.mat']);
    load([folder filesep 'timelock_oddball.mat']);
    
    a = standard.avg(effect_window_pos);
    b = oddball.avg(effect_window_pos);
    c = a-b;
    Pos.ERP_Diff_alltimechan(ii,:) =c;
    Pos.ERP_Diff(ii) = nanmean([a - b]);
    clear expected unexpected a b c

end

% Calculate Cohen's d
Pos.stdev_ERP_diff      = std(Pos.ERP_Diff);
Pos.mean_ERP_diff       = mean(Pos.ERP_Diff);
Pos.cohensd_ERP_diff    = Pos.mean_ERP_diff/Pos.stdev_ERP_diff;

% Then for the negative clustesr
effect_window_neg = stat_standard_oddball_clusstats.negclusterslabelmat==1;

% Calculate pairwise difference between conditions for each participant
for ii = 1:size(subjectlist_new,1)
    folder                  = [results filesep subjectlist_new{ii}];
    load([folder filesep 'timelock_standard.mat']);
    load([folder filesep 'timelock_oddball.mat']);
    
    a = standard.avg(effect_window_neg);
    b = oddball.avg(effect_window_neg);
    c = a-b;
    Neg.ERP_Diff_alltimechan(ii,:) =c;
    Neg.ERP_Diff(ii) = nanmean([a - b]);
    clear expected unexpected a b c

end

% Calculate Cohen's d
Neg.stdev_ERP_diff      = std(Neg.ERP_Diff);
Neg.mean_ERP_diff       = mean(Neg.ERP_Diff);
Neg.cohensd_ERP_diff    = Neg.mean_ERP_diff/Neg.stdev_ERP_diff;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Option 2: Determine at maximum effect size and at which channel/time it
% is maximal (upper bound)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determine maximum effect size and at which channel and time point Cohen's d is maximal
for t = 1:size(Pos.ERP_Diff_alltimechan,2)   
    Pos.cohensd_ERP_Diff_alltimechan(t) =(nanmean(Pos.ERP_Diff_alltimechan(:,t)))/(std(Pos.ERP_Diff_alltimechan(:,t)));
end

[Pos.cohensd_ERP_Diff_Max, Pos.idx]     = max(Pos.cohensd_ERP_Diff_alltimechan);
[Pos.row,Pos.col]                       = find(stat_standard_oddball_clusstats.posclusterslabelmat==1);
 
% Determine maximum effect size and at which channel and time point Cohen's d is maximal

for t = 1:size(Neg.ERP_Diff_alltimechan,2)   
    Neg.cohensd_ERP_Diff_alltimechan(t) =(nanmean(Neg.ERP_Diff_alltimechan(:,t)))/(std(Neg.ERP_Diff_alltimechan(:,t)));
end

[Neg.cohensd_ERP_Diff_Max, Neg.idx]     = min(Neg.cohensd_ERP_Diff_alltimechan);
[Neg.row,Neg.col]                       = find(stat_standard_oddball_clusstats.negclusterslabelmat==1);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Option 3: % Calculate effect size on rectangle around cluster results (lower bound)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate grand average, keeping information about individual
% participants
cfg                      = [];
cfg.channel              = 'all';
cfg.latency              = 'all';
cfg.parameter            = 'avg';
cfg.keepindividual       = 'yes';
grandavg_standard_all    = ft_timelockgrandaverage(cfg, standard_all{:});
grandavg_oddball_all     = ft_timelockgrandaverage(cfg, oddball_all{:});

% First for the largest positive cluster

% Determine time and channels of the largest positive cluster 
[Pos.row,Pos.col] = find(stat_standard_oddball_clusstats.posclusterslabelmat==1); % row = channel; col = time
idx_time_min = min(Pos.col);
idx_time_max = max(Pos.col);

% Estimate rectangular window around this cluster
rect_t_min = stat_standard_oddball_clusstats.time(Pos.col(idx_time_min));
rect_t_max = stat_standard_oddball_clusstats.time(Pos.col(idx_time_max));
rect_chan = stat_standard_oddball_clusstats.label(any(stat_standard_oddball_clusstats.mask(:,idx_time_min:idx_time_max),2));

% Calculate effect size (Cohen's d) withing this rectengular window
cfg                     = [];
cfg.channel             = rect_chan;
cfg.latency             = [rect_t_min rect_t_max];
cfg.avgoverchan         = 'yes';
cfg.avgovertime         = 'yes';
cfg.method              = 'analytic';
cfg.statistic           = 'cohensd';
cfg.ivar                = 1;
cfg.uvar                = 2;
Nsub                    = length(subjectlist_new);
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];

effect_rectangle_pos = ft_timelockstatistics(cfg, grandavg_standard_all, grandavg_oddball_all);
Pos.effect_rectangle = effect_rectangle_pos;

% Determine time and channels of the largest negative cluster 
[Neg.row,Neg.col] = find(stat_standard_oddball_clusstats.negclusterslabelmat==1);
idx_time_min = min(Neg.col);
idx_time_max = max(Neg.col);

% Estimate rectangular window around this cluster
rect_t_min = stat_standard_oddball_clusstats.time(Neg.col(idx_time_min));
rect_t_max = stat_standard_oddball_clusstats.time(Neg.col(idx_time_max));
rect_chan = stat_standard_oddball_clusstats.label(any(stat_standard_oddball_clusstats.mask(:,idx_time_min:idx_time_max),2));

% Calculate effect size (Cohen's d) withing this rectengular window
cfg                     = [];
cfg.channel             = rect_chan;
cfg.latency             = [rect_t_min rect_t_max];
cfg.avgoverchan         = 'yes';
cfg.avgovertime         = 'yes';
cfg.method              = 'analytic';
cfg.statistic           = 'cohensd';
cfg.ivar                = 1;
cfg.uvar                = 2;
Nsub                    = length(subjectlist_new);
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];

effect_rectangle_neg = ft_timelockstatistics(cfg, grandavg_standard_all, grandavg_oddball_all);
Neg.effect_rectangle = effect_rectangle_neg;

% Display results
fprintf('\n')
disp('~~~~~')
disp(['Contrast: Standard vs. Oddball: ']);
disp(['Effect size Cohens d of average positive cluster is ' num2str(Pos.cohensd_ERP_diff)])
disp(['Maximum effect size is ' num2str(Pos.cohensd_ERP_Diff_Max) ' at channel ' standard.label{Pos.row(Pos.idx)} ' and at time ' num2str(standard.time(Pos.col(Pos.idx))) ' sec']);
disp(['Effect size Cohens d of average on rectangle around positive cluster is ' num2str(effect_rectangle_pos.cohensd)])
disp('~~~~~')
fprintf('\n')

fprintf('\n')
disp('~~~~~')
disp(['Contrast: Standard vs. Oddball: '])
disp(['Effect size Cohens d of average negative cluster is ' num2str(Neg.cohensd_ERP_diff)])
disp(['Maximum effect size is ' num2str(Neg.cohensd_ERP_Diff_Max) ' at channel ' standard.label{Neg.row(Neg.idx)} ' and at time ' num2str(standard.time(Neg.col(Neg.idx))) ' sec']);
disp(['Effect size Cohens d of average on rectangle around negative cluster is ' num2str(effect_rectangle_neg.cohensd)])
disp('~~~~~')
fprintf('\n')

% Save the data
save(fullfile(output_dir, 'EffectSize.mat'), 'Neg', 'Pos');

%% Plot ERP timecourse of the channels with the maximal effect size

% Determine variability between participants
se_grandavg_standard = squeeze(nanstd(grandavg_standard_all.individual/sqrt(length(subjectlist_new))));
se_grandavg_oddball = squeeze(nanstd(grandavg_oddball_all.individual/sqrt(length(subjectlist_new))));

% Set color specifications
colour_code = {'b','r', 'k'};
shaded_area = {[0, 0, 1], [1 0 0], [0, 0, 0]};

% ERP of channel with maximum effect size of positive Cluster
figure;

% Condition 1
subplot(1,2,1)
plot(grandavg_standard.time,grandavg_standard.avg(Pos.row(Pos.idx),:),colour_code{1}, 'LineWidth', 1.5)
mean_standard = grandavg_standard.avg(Pos.row(Pos.idx),:);
se_standard = se_grandavg_standard(Pos.row(Pos.idx),:);
patch([grandavg_standard.time, fliplr(grandavg_standard.time)], [mean_standard-se_standard, fliplr(mean_standard+se_standard)],  shaded_area{1}, 'edgecolor', 'none', 'FaceAlpha', .3);

hold all

% Condition 2
plot(grandavg_oddball.time,grandavg_oddball.avg(Pos.row(Pos.idx),:),colour_code{2}, 'LineWidth', 1.5)
mean_oddball = grandavg_oddball.avg(Pos.row(Pos.idx),:);
se_oddball = se_grandavg_oddball(Pos.row(Pos.idx),:);
patch([grandavg_standard.time, fliplr(grandavg_standard.time)], [mean_oddball-se_oddball, fliplr(mean_oddball+se_oddball)],  shaded_area{2}, 'edgecolor', 'none', 'FaceAlpha', .3);

% Shaded area to indicate positive cluster
idx_pos_time = find(stat_standard_oddball_clusstats.posclusterslabelmat(Pos.row(Pos.idx),:)==1);
hold all
patch([grandavg_standard.time(idx_pos_time),fliplr(grandavg_standard.time(idx_pos_time))], [(ones(size(grandavg_standard.time(idx_pos_time),2),1)*-15)', fliplr((ones(size(grandavg_standard.time(idx_pos_time),2),1)*15)')], shaded_area{3}, 'edgecolor', 'none', 'FaceAlpha', .1)

xlabel('Time [s]');
ylabel('Amplitude [mV]');
ylim([-15 15])
line([-.5 1],[ 0 0], 'Color', [0 0 0],'LineStyle', ':')
title(['Maximum effect of positive cluster at channel ' grandavg_standard.label(Pos.row(Pos.idx))])

% ERP of channel with maximum effect size of negative Cluster
% Condition 1
subplot(1,2,2)
plot(grandavg_standard.time,grandavg_standard.avg(Neg.row(Neg.idx),:),colour_code{1}, 'LineWidth', 1.5)
mean_standard = grandavg_standard.avg(Neg.row(Neg.idx),:);
se_standard = se_grandavg_standard(Neg.row(Neg.idx),:);
patch([grandavg_standard.time, fliplr(grandavg_standard.time)], [mean_standard-se_standard, fliplr(mean_standard+se_standard)],  shaded_area{1}, 'edgecolor', 'none', 'FaceAlpha', .3);

hold all

% Condition 2
plot(grandavg_oddball.time,grandavg_oddball.avg(Neg.row(Neg.idx),:),colour_code{2}, 'LineWidth', 1.5)
mean_oddball = grandavg_oddball.avg(Neg.row(Neg.idx),:);
se_oddball = se_grandavg_oddball(Neg.row(Neg.idx),:);
patch([grandavg_standard.time, fliplr(grandavg_standard.time)], [mean_oddball-se_oddball, fliplr(mean_oddball+se_oddball)],  shaded_area{2}, 'edgecolor', 'none', 'FaceAlpha', .3);

% Shaded area to indicate negative cluster
idx_neg_time = find(stat_standard_oddball_clusstats.negclusterslabelmat(Neg.row(Neg.idx),:)==1);
hold all
patch([grandavg_standard.time(idx_neg_time),fliplr(grandavg_standard.time(idx_neg_time))], [(ones(size(grandavg_standard.time(idx_neg_time),2),1)*-15)', fliplr((ones(size(grandavg_standard.time(idx_neg_time),2),1)*15)')], shaded_area{3}, 'edgecolor', 'none', 'FaceAlpha', .1)

xlabel('Time [s]');
ylabel('Amplitude [mV]');
ylim([-15 15])
line([-.5 1],[ 0 0], 'Color', [0 0 0],'LineStyle', ':')
title(['Maximum effect of negative cluster at channel ' grandavg_standard.label(Neg.row(Neg.idx))])

%% Plot effect size topography highlighting cluster-based permutation test results

% Determine effect size for each channel x time pair
cfg                 = [];
cfg.parameter       = 'individual';
cfg.method          = 'analytic';
cfg.statistic       = 'cohensd';
cfg.ivar            = 1;
cfg.uvar            = 2;
num_sub             = length(standard_all);
cfg.design          = [ones(1,num_sub) ones(1,num_sub)*2 
    1:num_sub 1:num_sub];

effect_all_with_mask = ft_timelockstatistics(cfg, grandavg_standard_all, grandavg_oddball_all);

% Create mask to indicate clusters
effect_all_with_mask.mask = stat_standard_oddball_clusstats.mask;

cfg                 = [];
cfg.layout          = 'EEG1010.lay';
cfg.parameter       = 'cohensd';
cfg.maskparameter   = 'mask';
cfg.linecolor       = [0 0 0];
ft_multiplotER(cfg,effect_all_with_mask)

close all
