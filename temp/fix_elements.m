%%
load('lrb85315auto_annotation4.mat')
idxs = [];
for i = 1:numel(keys)
tokens = regexp(keys{i},'_','split');
idxs = [idxs; str2num(tokens{2})];
end
loc = find(idxs == 4876)
%%
elements{loc}
%%
locs = find(elements{loc}.segType == 103)
%%
elements{loc}.segType(locs) = 405; 