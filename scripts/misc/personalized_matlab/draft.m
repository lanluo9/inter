tt = fieldnames(input);
index = cellfun(@(x) any(contains(x, 'On')),tt); sum(index)
id = find(index > 0);
for i = 1 : length(id)
    fprintf(['input.', tt{id(i)}])
    fprintf('\n')
end

% input.cStimOneOff