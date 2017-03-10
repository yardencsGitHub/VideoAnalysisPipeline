function [unique_words, locations] = words_in_batch(raster)
% Takes a binary raster [bits x words] and returns the decimal rep. of the
% unique words in the raster and a cell array (locations) of their
% locations in the raster.
n = size(raster,1);
words = raster'*(2.^[(n-1):-1:0]');
unique_words = unique(words);
locations = {};
parfor word_cnt = 1:numel(unique_words)
    locations{word_cnt} = find(words == unique_words(word_cnt));
end
