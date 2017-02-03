% Take the movie Y [width x height x frames] and align frames to min
% consider changing to imregister
function [Ytmp,d]  = AlignSingleTrial_sequential(Y)
[optimizer,metric] = imregconfig('monomodal');
Fr = double(Y);
%Fr = Fr - imgaussfilt(Fr,15);
Ytmp = double(Y);
d={}
Npixels = numel(Y(:,:,1));
for frm = size(Y,3):-1:2
    display(frm)
    X = Fr(:,:,frm);
    tform = imregtform(Fr(:,:,frm-1), X,'rigid',optimizer,metric);
    d = {tform d{:}};
    MOVING_REG = imwarp(Fr(:,:,frm-1),tform,'OutputView',imref2d(size(Y(:,:,1))));
    %MOVING_REG = imregister(Fr(:,:,frm-1), X, 'rigid', optimizer,metric);
    Fr(:,:,frm-1) = MOVING_REG;
    Ytmp(:,:,frm-1) = imwarp(Ytmp(:,:,frm-1),tform,'OutputView',imref2d(size(Y(:,:,1))));
end
    
  