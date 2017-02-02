% Take the movie Y [width x height x frames] and align frames to min
% consider changing to imregister
function [Ytmp,d]  = AlignSingleTrial(Y)
%X = (mean(double(Y(:,:,20:end)),3)); 
%X = (min(double(Y(:,:,20:end)),[],3)); %25:200,50:225
%X = (quantile(double(Y(:,:,20:end)),0.25,3));
%X = imgaussfilt(X,10);
X =  double(Y(:,:,11));
%X = X-min(X(:));
%X = X/max(X(:));

[optimizer,metric] = imregconfig('monomodal');
Fr = double(Y);
%Fr = Fr - min(Fr(:));
%Fr = Fr / max(Fr(:));
Ytmp = Y;
d=zeros(size(Y,3),2);
Npixels = numel(Y(:,:,1));
parfor frm = 1:size(Y,3)
    %disp(frm)
    %Fr = double(squeeze(Y(:,:,frm)));
    %Fr(:,:,frm) = Fr(:,:,frm) - min(reshape(Fr(:,:,frm),1,Npixels));
    %Fr(:,:,frm) = Fr(:,:,frm) / max(reshape(Fr(:,:,frm),1,Npixels));
    MOVING_REG = imregister(Fr(:,:,frm), X, 'rigid', optimizer,metric);
    %tform = imregtform(Fr(:,:,frm),X,'rigid',optimizer,metric);
%     d(frm,:) = tform.T(3,1:2);
%     if (max(abs(d(frm,:)))) < 20
%         Ytmp(:,:,frm) = imwarp(Y(:,:,frm),tform,'OutputView',imref2d(size(squeeze(Y(:,:,1)))));
%     else
%         Ytmp(:,:,frm) = Y(:,:,frm); 
%         display(frm);
%     end
    Ytmp(:,:,frm) = MOVING_REG;
    
end
    
    