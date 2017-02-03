function r = sbxalign2(Y,idx)
[optimizer,metric] = imregconfig('monomodal');
if(length(idx)==1)
     A = Y(:,:,idx);
     r.m = A; 
     r.T = {eye(3)};
     r.idx = idx;
else
 idx0 = idx(1:floor(end/2)); % split into two groups
 idx1 = idx(floor(end/2)+1 : end);
 
 r0 = sbxalign2(Y,idx0); % align each group
 r1 = sbxalign2(Y,idx1);
 tform = imregtform(r0.m, r1.m,'rigid',optimizer,metric);
 shift = sqrt(sum(tform.T(3,1:2).^2));
 angle = atan(abs(tform.T(1,2)/tform.T(1,1)));
 %if ((shift > 10) || angle/pi*180 > 10)
  %   tform = affine2d(eye(3));
 %end
 MOVING_REG = imwarp(r0.m,tform,'OutputView',imref2d(size(r0.m)));
 r.m = (MOVING_REG+r1.m)/2;
 if iscell(r0.T)
    q = cellfun(@(x)(tform.T*x),r0.T,'UniformOutput',false);
    if iscell(r1.T)
        r.T = {q{:} r1.T{:}};
    else
        r.T = {q{:} r1.T};
    end
 else
    q = tform.T*r0.T;
    if iscell(r1.T)
        r.T = {q r1.T{:}};
    else
        r.T = {q r1.T};
    end
 end
 r.idx = [r0.idx;r1.idx];
end