function write_matrix_to_file(data, filename)
data = full(data); % ensure that data is not sparse format
fid = fopen(filename, 'w');
rowN = size(data, 1);
for i = 1: rowN
    fprintf(fid, '%.17g	', data(i, :));
    fprintf(fid, '\n');
end
fclose(fid);