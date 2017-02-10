function movie = Mat2Mov(inmat)
% used for a quick checking of a 3d matrix
Xt = permute(inmat,[1 2 4 3]); % 4D matrix
movie = immovie(uint8(Xt),gray(256)); 
%implay(movie);