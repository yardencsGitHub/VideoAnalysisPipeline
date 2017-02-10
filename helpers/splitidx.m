function res = splitidx(b,numparts)
l = length(b);
d = floor(l/numparts);
r = mod(l,d);
a = 1:l;
borders = [a(1:d:end-d+1)' a(d:d:end)'];
if (r>0)
    borders = [borders; borders(end)+1 borders(end)+r];
end
res = cell(numparts,1);
for i = 1:numparts
    res{i} = b(borders(i,1):borders(i,2));
end