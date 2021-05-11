function trl = trialfun_BeeG(cfg)

% Read relevant info
event = ft_read_event(cfg.dataset);
hdr   = ft_read_header(cfg.dataset);

% Initiate the output struct
counttrial = 0;

% Set up the cell arrays contaning marker values of each condition
expected_array = {'S  1','S  2','S  3','S  4','S  5','S  6','S  7','S  8'};
post_update_array = {'S 11','S 12','S 13','S 14','S 15','S 16','S 17','S 18'};
post_no_update_array = {'S 21','S 22','S 23','S 24','S 25','S 26','S 27','S 28'};
locations = [0, 45, 90, 135, 180, 225, 270, 315]';

% Now let's loop over all events
for ii = 1:size(event,2)
  % First we test if it is an event of interest
  if strcmp(event(ii).type, cfg.trialdef.eventtype)
    % Then it is a stimulus of interest
    
    % Test which stimulus it is
    if strcmp(event(ii).value, 'S  9')
      % This is a blank screen displayed before the fixation cross
      counttrial = counttrial + 1;
      begsample(counttrial, :) = event(ii).sample; % Start sample number of the event
      endsample(counttrial, :) = event(ii).sample + round(0.5 * hdr.Fs);  % the stimulus is 500 ms long
      offset(counttrial, :) = 0;
      marker(counttrial, :) = {event(ii).value};
      stimulus(counttrial, : ) = {'blank screen'};
      location_bee(counttrial, :) = nan;
    elseif strcmp(event(ii).value, 'S 10')
      % This is a fixation cross
      counttrial = counttrial + 1;
      begsample(counttrial, :) = event(ii).sample; % Start sample number of the event
      endsample(counttrial, :) = event(ii).sample + round(1 * hdr.Fs);  % the stimulus is 1000 ms
      offset(counttrial, :) = 0;
      marker(counttrial, :) = {event(ii).value};
      stimulus(counttrial, : ) = {'fixation cross'};
      location_bee(counttrial, :) = nan;
    elseif any(strcmp(event(ii).value, expected_array))
      % Then it is an expected bee
      counttrial = counttrial + 1;
      begsample(counttrial, :) = event(ii).sample; % Start sample number of the event
      endsample(counttrial, :) = event(ii).sample + round(1.5 * hdr.Fs);  % the stimulus is 1000 ms
      offset(counttrial, :) = 0;
      marker(counttrial, :) = {event(ii).value};
      stimulus(counttrial, : ) = {'bee'};
      index = find(strcmp(event(ii).value, expected_array));
      location_bee(counttrial, :) = locations(index);
    elseif  any(strcmp(event(ii).value, post_update_array))
      % Then it is a bee presented after an update cue
      counttrial = counttrial + 1;
      begsample(counttrial, :) = event(ii).sample; % Start sample number of the event
      endsample(counttrial, :) = event(ii).sample + round(1.5 * hdr.Fs);  % the stimulus is 1000 ms
      offset(counttrial, :) = 0;
      marker(counttrial, :) = {event(ii).value};
      stimulus(counttrial, : ) = {'post update-cue bee'};
      index = find(strcmp(event(ii).value, post_update_array));
      location_bee(counttrial, :) = locations(index);
    elseif any(strcmp(event(ii).value, post_no_update_array))
      % Then it is a bee presented after a no-update cue
      counttrial = counttrial + 1;
      begsample(counttrial, :) = event(ii).sample; % Start sample number of the event
      endsample(counttrial, :) = event(ii).sample + round(1.5 * hdr.Fs);  % the stimulus is 1000 ms
      offset(counttrial, :) = 0;
      marker(counttrial, :) = {event(ii).value};
      stimulus(counttrial, : ) = {'post no-update-cue bee'};
      index = find(strcmp(event(ii).value, post_no_update_array));
      location_bee(counttrial, :) = locations(index);
    elseif strcmp(event(ii).value, 'S 19')
      % Then it is an update cue
      counttrial = counttrial + 1;
      begsample(counttrial, :) = event(ii).sample; % Start sample number of the event
      endsample(counttrial, :) = event(ii).sample + round(1.5 * hdr.Fs);  % the stimulus is 1000 ms
      offset(counttrial, :) = 0;
      marker(counttrial, :) = {event(ii).value};
      stimulus(counttrial, : ) = {'update-cue'};
      location_bee(counttrial, :) = nan;
    elseif strcmp(event(ii).value, 'S 20')
      % Then it is a no-update cue
      counttrial = counttrial + 1;
      begsample(counttrial, :) = event(ii).sample; % Start sample number of the event
      endsample(counttrial, :) = event(ii).sample + round(1.5 * hdr.Fs);  % the stimulus is 1000 ms
      offset(counttrial, :) = 0;
      marker(counttrial, :) = {event(ii).value};
      stimulus(counttrial, : ) = {'no-update-cue'};
      location_bee(counttrial, :) = nan;
    end
  end
end

trl = table(begsample, endsample, offset, marker, stimulus, location_bee);