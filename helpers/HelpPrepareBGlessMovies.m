%% To prepare BG subtracted movies ...
b1 = quantile(Y(:,:,n_del_frames+1:end),0.1,3);
b1 = (b1-min(b1(:)))/max(b1(:));
b1 = b1/median(b1(:));


bt = median(reshape(Y,rows*columns,frames));
b2 = bsxfun(@times,b1,reshape(bt,1,1,numel(bt)));
c = Y-b2; c = c - min(min(min(c(:,:,n_del_frames+1:end))));
%%
c1 = c(:,:,5:end);
c1 = c1/max(c1(:))*255;

%%
dff_temp = bsxfun(@rdivide,bsxfun(@minus,c1,quantile(c1,0.05,3)),quantile(c1,0.05,3));