%%
if ~exist('MaxProj','dir')
    mkdir('MaxProj');
end
map = [0:0.01:1]'*ones(1,3);
filt_rad = 50; filt_sigma = 45; n_del_frames = 0;
h = fspecial('gaussian',filt_rad,filt_sigma);
h1 = fspecial('disk',2); 
% OutputDIR = '/Users/yardenc/Documents/Experiments/Imaging/CanaryData/lrb853_15/movs/MaxProj';
% DIR = '/Volumes/home/Data/Imaging/lrb853_15/RawData';
FILES = dir(fullfile(DIR,'*.mat'));
ord = [];
for i = 1:numel(FILES)
    tokens = regexp(FILES(i).name,'_','split');
    ord = [ord; str2num(tokens{3})];
end
[locs,indx] = sort(ord);
startloc = min(find(locs >= startfrom));
FILES = FILES(indx);
for fnum = startloc:numel(FILES)
    load(fullfile(DIR,FILES(fnum).name));
    tokens = regexp(FILES(fnum).name,'_','split');
    %Y = Y(:,:,n_del_frames+1:end);
    Y = imresize(vidMat(:,:,n_del_frames+1:end),0.5);
    Y = Y - min(Y(:));
    Y = Y / max(Y(:));
    Y = convn(Y,reshape([1 1 1 1 1] / 5, 1, 1, []), 'same');
    base = imfilter(Y,h,'circular','replicate');
    Y = (Y - base);
    Y = bsxfun(@minus,Y,quantile(Y,0.05,3)); %bsxfun(@rdivide,,quantile(Y,0.05,3))
    maxY = max(Y,[],3);
    maxY(maxY(:) < 0.05) = 0;
    out_filename = char(join(tokens(2:end),'_'));
    out_filename(end-2:end)='png';
    %[FILES(fnum).name(9:end-3) 'png'];
    imwrite(uint8(256*maxY/max(max(maxY(:)),0.4)),gray(256),fullfile(OutputDIR,out_filename));
    display([fnum numel(FILES)]);
end
%% HISTORICAL optional: remove initial filename segment
% cd(OutputDIR);
% for i = 1:numel(d)
%     movefile(d(i).name, d(i).name(9:end));
% end
