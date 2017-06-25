cd('/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/Images');
old = csvread('Discarded_or_noise.csv');
cd('/Users/yardenc/Documents/Experiments/Imaging/Data/CanaryData/lrb853_15/movs/wav/annotated/images');
d = dir('*.png');
nums = zeros(numel(d),1);
new = old;
for fnum = 1:numel(d)
    tokens = regexp(d(fnum).name,'_','split');
    num = str2num(tokens{2});
    if ~ismember(num,old)
        I = imread(d(fnum).name);
        h = figure;
        imagesc(I);
        res = input(['File # ' tokens{2} ' noise?:']);
        if (res == 1)
            new = [new; num];
            display('NOISE!');
        end
        hgclose(h);
    end
end
