function write_json(filename, json)

ft_info('writing ''%s''\n', filename);
json = remove_empty(json);
% json = sort_fields(json) % Let's leave this out so that the order will remain as we indicate it here in the code, not alphabetically ordered;
json = ft_struct2char(json); % convert strings into char-arrays
ft_hastoolbox('jsonlab', 1);
% see also the output_compatible helper function
% write nan as 'n/a'
% write boolean as True/False
str = savejson('', json, 'NaN', '"n/a"', 'ParseLogical', true);
% fid = fopen_or_error(filename, 'w');
fid = fopen(filename, 'w');
fwrite(fid, str);
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUBFUNCTION
function s = remove_empty(s)
if isempty(s)
  return
elseif isstruct(s)
  fn = fieldnames(s);
  fn = fn(structfun(@isempty, s));
  s = removefields(s, fn);
elseif istable(s)
  remove = false(1,size(s,2));
  for i=1:size(s,2)
    % find columns that are non-numeric and where all elements are []
    remove(i) = ~isnumeric(s{:,i}) && all(cellfun(@isempty, s{:,i}));
  end
  s = s(:,~remove);
end
