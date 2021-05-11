%% Define neighbours:
%
% These scripts and the data in BIDS format are part of Meyer, M., Lamers, D.,
% Kayhan, E., Hunnius, S., & Oostenveld, R. (2021) Fostering reproducibility in
% developmental EEG research by using BIDS, cluster-based permutation tests and
% reporting effectsizes (in preparation)
%
% The infant EEG dataset is originally described in Kayhan, E., Meyer, M., O'Reilly,
% J. X., Hunnius, S., & Bekkering, H. (2019). Nine-month-old infants update their
% predictive models of a changing environment. Developmental cognitive neuroscience,
% 38, 100680.)
%
% The neighbours are used to interpolate bad channels and for cluster-based
% statistic. The neighbours are specific for the EEG setup and electrode layout used
% in this study, which comprises 32 channels that were placed according to the 10-20
% standard.
%
% The code in this script is not meant to be executed as-is, but is provided for
% reference. We selected the neighbours in an interactive way in the 1st part. One
% pair was missed, which is dealt with in the 2nd refinement part. Also, neighbours
% are made symmetric in the 2nd part.

if isfile('selected_neighbours.mat')
  load(fullfile(scripts, 'selected_neighbours.mat'));
  do_initial    = false;
  do_refinement = false;
else
  do_initial    = true;
  do_refinement = true;
end

%% Read the subject info from the BIDS dataset

t           = readtable([bidsroot filesep 'participants.tsv'], 'FileType', 'text');
subjectlist = t.participant_id;

%% Read the channels used in the EEG study from the first (representative) subject

sub     = subjectlist{1};
dataset = [bidsroot filesep sub filesep 'eeg' filesep sub '_task-audiovisual_eeg.vhdr' ];
hdr     = ft_read_header(dataset);
label   = hdr.label;

% Read template 3D sensor positions for the standard 1020 electrode placement system
elec = ft_read_sens('standard_1020.elc', 'senstype', 'eeg');

% Now loop through all used channels, to create a new layout struct
% containing only those channels used in this particular study

for ii = 1:length(label)
  index                         = strcmp(label{ii}, elec.label);
  selected_elec.chanpos(ii, :)  = elec.chanpos(index, :);
  selected_elec.chantype{ii, 1} = elec.chantype{index};
  selected_elec.chanunit{ii, 1} = elec.chanunit{index};
  selected_elec.elecpos(ii, :)  = elec.elecpos(index, :);
  selected_elec.label{ii, 1}    = elec.label{index};
  selected_elec.type            = elec.type;
  selected_elec.unit            = elec.unit;
end

%% Draw lines in the figure to connect channels and create a neighbours struct

if do_initial
  
  % Plot the selected electrodes in 3D
  fig = ft_plot_sens(selected_elec, 'label', 'label');
  
  % Loop through all channels
  for cc = 1:length(selected_elec.label)
    
    % Find the corresponding label
    channel_to_connect = selected_elec.label{cc};
    
    % Allow the user to provide a cell array with the channels to which it is connected
    text = ['give a cell-array with the neighbours to ' channel_to_connect ' : '];
    connected_channels = input(text);
    
    % Create a line between the selected channel and the channels to which it is connected
    connect1 = find(strcmp(selected_elec.label, channel_to_connect));
    for ci = 1:length(connected_channels)
      connect2    = find(strcmp(selected_elec.label, connected_channels(ci)));
      linepoints  = [selected_elec.chanpos(connect1, :); selected_elec.chanpos(connect2, :)];
      line(linepoints(:,1), linepoints(:,2), linepoints(:,3));
    end
    
    % Create the neighbours struct
    selected_neighbours(cc).label       = channel_to_connect;
    selected_neighbours(cc).neighblabel = connected_channels;
  end % for each label
  
  % Save the neighbours struct
  save(fullfile(scripts, 'selected_neighbours.mat'), 'selected_neighbours');
  savefig(fig, fullfile(scripts, 'selected_neighbours_plot'));
  
end % if do initial

% When the above part of the code is finished, comment it out to continue with some
% manual optimization in the next step

%%  Make some manual changes if necessary

if isfile('selected_neighbours.mat')
  load('selected_neighbours.mat');
end

if isfile('selected_neighbours_plot.fig')
  openfig('selected_neighbours_plot.fig');
end

if do_refinement
  % Add a line between two points that were missed by the previous step:
  connect1    = find(strcmp(selected_elec.label, 'Fp1'));
  connect2    = find(strcmp(selected_elec.label, 'Fp2'));
  linepoints  = [selected_elec.chanpos(connect1, :); selected_elec.chanpos(connect2, :)];
  line(linepoints(:,1), linepoints(:,2), linepoints(:,3));
  
  % And add the connection to the neighours struct
  selected_neighbours(1).neighblabel  = [selected_neighbours(1).neighblabel, {'Fp2'}];
  selected_neighbours(2).neighblabel  = [selected_neighbours(2).neighblabel, {'Fp1'}];
  
  % Then do a sanity check, make sure that neighbours are fully bidirectional
  % I.e. if Cz has FCz as a neighbour, then FCz should have Cz as a neighbour
  
  for ii = 1:length(selected_neighbours)
    channel = selected_neighbours(ii).label;
    for tt  = 1:length(selected_neighbours(ii).neighblabel)
      idx = find(strcmp(label, selected_neighbours(ii).neighblabel{tt}));
      if sum(strcmp(selected_neighbours(idx).neighblabel, channel))==0
        selected_neighbours(idx).neighblabel = [selected_neighbours(idx).neighblabel, {channel}];
      end
    end
  end
  
  % Finally, save this plot and new neighbours struct
  save(fullfile(scripts, 'selected_neighbours.mat'), 'selected_neighbours');
  savefig(fig, fullfile(scripts, 'selected_neighbours_plot'));
  
end % if do refinement

%% And replot so that all lines are now corrected

close all

fig = ft_plot_sens(selected_elec, 'label', 'label');

for ii = 1:length(selected_neighbours)
  channel         = selected_neighbours(ii).label;
  connect1         = find(strcmp(selected_elec.label, channel));
  for tt          = 1:length(selected_neighbours(ii).neighblabel)
    connect2    = find(strcmp(selected_elec.label, selected_neighbours(ii).neighblabel{tt}));
    linepoints  = [selected_elec.chanpos(connect1, :); selected_elec.chanpos(connect2, :)];
    line(linepoints(:,1), linepoints(:,2), linepoints(:,3));
  end
end
