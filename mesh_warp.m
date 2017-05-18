function [out_x, out_y, Hg_mask, Hmdlt] = mesh_warp(grid_size, imgw, imgh, warp_type, T, i, ref, tar, H, c1para, mdltpara, uvmesh_ON)

% mesh by meshgrid
x_list = linspace(1, imgw, ceil(imgw/grid_size));
y_list = linspace(1, imgh, ceil(imgh/grid_size));
[X0, Y0] = meshgrid(x_list, y_list);
[X, Y] = apply_transform(X0, Y0, T{i});

% uv mesh
if ~exist('uvmesh_ON', 'var'), uvmesh_ON = 0; end
if uvmesh_ON
    theta = c1para.theta;
    R_uv_to_xy = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];
    [X, Y] = apply_transform(X, Y, R_uv_to_xy);
end

if strcmp(warp_type, 'ours')
    [tmpx, tmpy] = apply_transform(X, Y, H{ref, i});
    [out_x, out_y] = c1_warp_ver2(tmpx, tmpy, c1para);
elseif strcmp(warp_type, 'hom')
    [out_x, out_y] = apply_transform(X, Y, H{tar, i});
else
    fprintf('frank: wrong warp type\n'); pause
end