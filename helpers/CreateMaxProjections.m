%%
map = [0:0.01:1]'*ones(1,3);
filt_rad = 50; filt_sigma = 45; n_del_frames = 5;
h = fspecial('gaussian',filt_rad,filt_sigma);
h1 = fspecial('disk',2); 
OutputDIR = '/Users/yardenc/Documents/Experiments/Imaging/CanaryData/lrb853_15/movs/MaxProj';
DIR = '/Volumes/home/Data/Imaging/lrb853_15/RawData';
FILES = dir(fullfile(DIR,'*.mat'));
for fnum = 1:numel(FILES)
    load(fullfile(DIR,FILES(fnum).name));
    Y = Y(:,:,n_del_frames+1:end);
    Y = imresize(vidMat,0.5);
    Y = Y - min(Y(:));
    Y = Y / max(Y(:));
    Y = convn(Y,reshape([1 1 1 1 1] / 5, 1, 1, []), 'same');
    base = imfilter(Y,h,'circular','replicate');
    Y = (Y - base);
    Y = bsxfun(@minus,Y,quantile(Y,0.05,3)); %bsxfun(@rdivide,,quantile(Y,0.05,3))
    maxY = max(Y,[],3);
    maxY(maxY(:) < 0.05) = 0;
    
    imwrite(uint8(256*maxY/max(max(maxY(:)),0.4)),gray(256),fullfile(OutputDIR,[FILES(fnum).name(1:end-3) 'png']));
    display([fnum numel(FILES)]);
end
%% optional: remove initial filename segment
cd(OutputDIR);
for i = 1:numel(d)
    movefile(d(i).name, d(i).name(9:end));
end
