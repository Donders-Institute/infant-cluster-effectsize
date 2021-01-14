function json = read_json(filename)
ft_info('reading ''%s''\n', filename);
ft_hastoolbox('jsonlab', 1);
json = loadjson(filename);
json = ft_struct2char(json); % convert strings into char-arrays