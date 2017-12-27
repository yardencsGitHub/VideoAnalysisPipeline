%%
Y = vidMat;
%% To prepare BG subtracted movies ...
% The movie Y is a 3d matrix rows x columns x frames

% The movie Y isa 3d matrix rows x columns x frames

% code doesn't deal with dropped frames and \ or motion correction (frame
% registration)
n_del_frames = 6; % How many frames to disregard at the video's onset
[rows,columns,frames] = size(Y);
%%
% calculate a single background frame - b1
b1 = quantile(Y(:,:,n_del_frames+1:end),0.1,3);
b1 = (b1-min(b1(:)))/(max(b1(:))-min(b1(:)));
b1 = b1/median(b1(:));

% calculate a backroung time series as the median of each frame
bt = median(reshape(Y,rows*columns,frames));
% create a BG movie by scaling the BG frame (b1) by the BG time series (bt) and get b2, a very meaningful name. 
b2 = bsxfun(@times,b1,reshape(bt,1,1,numel(bt)));
% subtract the BG movie, b2, from the original, Y, and get c.
c = Y-b2; 

%% remove a smoothed base
filt_rad = 150; filt_sigma = 145; % highpass filter 
h = fspecial('gaussian',filt_rad,filt_sigma);
base = imfilter(Y,h,'circular','replicate');
c = Y-base;
%c = c-min(c(:));
% the next steps are not needed unless planning to create a movie or a df/f
%% make sure all pixels after fram n_del_frames are positive
c = c - min(min(min(c(:,:,n_del_frames+1:end))));
%% optional, crop first 'n_del_frames' frames and normalize
c1 = c(:,:,n_del_frames:end);
%c1 = c1/max(c1(:))*255;

%% calculate df/f using each pixel's 5% quantile as baseline
dff_temp = bsxfun(@rdivide,bsxfun(@minus,c1,quantile(c1,0.05,3)),quantile(c1,0.05,3));
c = Y-b2; c = c - min(min(min(c(:,:,n_del_frames+1:end))));
%%
%c1 = c(:,:,5:end);
c1 = c/max(c(:))*255;

%%
dff_temp = bsxfun(@rdivide,bsxfun(@minus,c1,quantile(c1,0.05,3)),quantile(c1,0.05,3));
%%
dff = dff_temp*150;
dff(dff(:)<0) = 0;
dff(dff(:)>255) = 255;
write_mat_2_moviefile(dff,['Example_8076' '.mp4'],'MPEG-4',30);

