%% go over a list of images and mark the noise ones
DIR = pwd;
startfrom = 6427;

FILES = dir('*.gif');
ord = [];
for i = 1:numel(FILES)
    tokens = regexp(FILES(i).name,'_','split');
    ord = [ord; str2num(tokens{2})];
end
[locs,indx] = sort(ord);
startloc = find(locs == startfrom);
FILES = FILES(indx);
old = [];
noisefiles = [];
songfiles = [];
for fnum = startloc:numel(FILES)
    tokens = regexp(FILES(fnum).name,'_','split');
    num = str2num(tokens{2});
    if ~ismember(num,old)
        I = imread(FILES(fnum).name);
        h = figure;
        imagesc(I);
        set(h,'Position',[ 24         278        1256         420]);
        res = input(['File # ' tokens{2} ' noise?:']);
        if (res == 1)
            noisefiles = [noisefiles; num];
            display('NOISE!');
        else
            songfiles = [songfiles; num];
        end
        hgclose(h);
    end
end

%%
keys = {};
for fnum = startloc:numel(FILES)
    tokens = regexp(FILES(fnum).name,'_','split');
    num = str2num(tokens{2});
    if ismember(num,songfiles)
        keys = {keys{:} [FILES(fnum).name(1:end-3) 'mat']};
    end
end


%% historical
% cd('/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/Images');
% old = csvread('Discarded_or_noise.csv');
% cd('/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/movs/wav/annotated/images');
% d = dir('*.png');