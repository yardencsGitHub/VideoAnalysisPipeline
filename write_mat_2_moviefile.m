function write_mat_2_moviefile(mtx,movfile,movtype,fs)
movie = Mat2Mov(mtx);
v = VideoWriter(movfile,movtype);
open(v)
writeVideo(v,movie);
close(v)
