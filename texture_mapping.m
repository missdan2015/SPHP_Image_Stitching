function [c1out, c1omask, T] = texture_mapping(in_name, imgw, imgh, img_grid_size, T, warp_type, H, ref, tar, c1para, mdltpara, drawmesh, bgcolor)

mapping_location = 'texture_mapping/';

img_n = numel(in_name);
c1out = cell(img_n, 1);
c1omask = cell(img_n, 1);
inmesh_X0 = cell(img_n, 1);     inmesh_Y0 = cell(img_n, 1);
outmesh_X = cell(img_n, 1);     outmesh_Y = cell(img_n, 1);
curr_minx = zeros(img_n, 1);    curr_miny = zeros(img_n, 1);
curr_maxx = zeros(img_n, 1);    curr_maxy = zeros(img_n, 1);

for i = 1 : img_n
    % for texture mapping
    meshW = ceil(imgw(i) / img_grid_size);
    meshH = ceil(imgh(i) / img_grid_size);
    X0list = linspace(1, imgw(i), meshW);
    Y0list = linspace(1, imgh(i), meshH);
    [inmesh_X0{i}, inmesh_Y0{i}] = meshgrid(X0list, Y0list);
    
    % compute outmesh
    [outmesh_X{i}, outmesh_Y{i}] = mesh_warp(img_grid_size, imgw(i), imgh(i), warp_type, T, i, ref, tar, H, c1para, mdltpara);
    
    curr_minx(i) = min(outmesh_X{i}(:)); curr_miny(i) = min(outmesh_Y{i}(:));
    curr_maxx(i) = max(outmesh_X{i}(:)); curr_maxy(i) = max(outmesh_Y{i}(:));
end % up to now, the coord. system is x towards right, y towards up
all_minx = min(curr_minx); all_miny = min(curr_miny);
all_maxx = max(curr_maxx); all_maxy = max(curr_maxy);
% fprintf('x from %.2f ~ %.2f\n', all_minx, all_maxx);
% fprintf('y from %.2f ~ %.2f\n', all_miny, all_maxy);

%% reorganize
% change of coord. s.t. (1) the upperleft of bounding box of composite corresponds to (0,0)
%                       (2) x towards right, y-axis towards down
T1 = [1, 0, -all_minx; 0, -1, all_maxy; 0, 0, 1];
for i = 1 : img_n
    [outmesh_X{i}, outmesh_Y{i}] = apply_transform(outmesh_X{i}, outmesh_Y{i}, T1);
    % recompute min & max value
    curr_minx(i) = min(outmesh_X{i}(:)); curr_miny(i) = min(outmesh_Y{i}(:));
    curr_maxx(i) = max(outmesh_X{i}(:)); curr_maxy(i) = max(outmesh_Y{i}(:));
end % now the coord. system is x towards right, y towards down. bounding box aligns with both axis
all_minx = min(curr_minx); all_miny = min(curr_miny);
all_maxx = max(curr_maxx); all_maxy = max(curr_maxy);
bbarea = (all_maxx-all_minx)*(all_maxy-all_miny); % the area of bounding box
%if all_minx~=0 || all_miny~=0, fprintf('%.12f %.12f\n', all_minx, all_miny); pause; end

% scale down to avoid memory issue
max_size = 600*800;
scl = sqrt(max_size/bbarea);
T2 = [scl, 0, 0; 0, scl, 0; 0, 0, 1];
for i = 1 : img_n
    outmesh_X{i} = outmesh_X{i} * scl;
    outmesh_Y{i} = outmesh_Y{i} * scl;
    % recompute min & max value
    curr_minx(i) = min(outmesh_X{i}(:)); curr_miny(i) = min(outmesh_Y{i}(:));
    curr_maxx(i) = max(outmesh_X{i}(:)); curr_maxy(i) = max(outmesh_Y{i}(:));
end
all_minx = min(curr_minx); all_miny = min(curr_miny);
all_maxx = max(curr_maxx); all_maxy = max(curr_maxy);

% numerical post-processing (set near zero to be zero)
for i = 1 : img_n
    outmesh_X{i}( outmesh_X{i}<=(1e-7) ) = 0;
    outmesh_Y{i}( outmesh_Y{i}<=(1e-7) ) = 0;
end

% adding borders by translating the composite
tr_x = 5;
tr_y = 5;
T3 = [1, 0, tr_x; 0, 1, tr_y; 0, 0, 1];
for i = 1 : img_n
    outmesh_X{i} = outmesh_X{i} + tr_x;
    outmesh_Y{i} = outmesh_Y{i} + tr_y;
end

outimg_w = ceil(all_maxx - all_minx) + 1 + 2*tr_x;
outimg_h = ceil(all_maxy - all_miny) + 1 + 2*tr_y;


%% opencv
for i = 1 : img_n
    write_matrix_to_file([inmesh_X0{i}; inmesh_Y0{i}], [mapping_location 'inmesh.txt']);
    write_matrix_to_file([outmesh_X{i}; outmesh_Y{i}], [mapping_location 'outmesh.txt']);

    eval(['cd ' mapping_location])
    dos(sprintf('texture_mapping.exe ../%s %d %d %d %d %s %s %s %d %s', ...
        in_name{i}, meshW, meshH, outimg_w, outimg_h, 'inmesh.txt', 'outmesh.txt', 'out.jpg', drawmesh, bgcolor));
    c1out{i} = im2double(imread('output.jpg'));
    c1omask{i} = im2double(imread('out_mask.png'));
    cd ../
end

T4 = [1, 0, 1; 0, 1, 1; 0, 0, 1]; % translation by 1 pixels because matlab coordinate starts from 1

T = T4*T3*T2*T1;