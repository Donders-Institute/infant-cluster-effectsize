function write_tsv(filename, tsv)
ft_info('writing ''%s''\n', filename);
fn = tsv.Properties.VariableNames;
for i=1:numel(fn)
  % write [] as 'n/a'
  % write nan as 'n/a'
  % write boolean as 'True' or 'False'
  tsv.(fn{i}) = output_compatible(tsv.(fn{i}));
end
writetable(tsv, filename, 'Delimiter', 'tab', 'FileType', 'text');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUBFUNCTION
function val = output_compatible(val)
if istable(val)
  fn = val.Properties.VariableNames;
  for i=1:numel(fn)
    val.(fn{i}) = output_compatible(val.(fn{i}));
  end
elseif iscell(val)
  % use recursion to make all elements compatible
  val = cellfun(@output_compatible, val, 'UniformOutput', false);
elseif isnumeric(val) && numel(val)>1 && any(isnan(val))
  % convert and use recursion to make all elements compatible
  val = num2cell(val);
  val = cellfun(@output_compatible, val, 'UniformOutput', false);
else
  % write [] as 'n/a'
  % write nan as 'n/a'
  % write boolean as 'True' or 'False'
  if isempty(val)
    val = 'n/a';
  elseif isnan(val)
    val = 'n/a';
  elseif islogical(val)
    if val
      val = 'True';
    else
      val = 'False';
    end
  end
end