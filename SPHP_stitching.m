function SPHP_stitching(in_name, edge_list, ref, tar, warp_type, zeroR_ON)

img_n = size(in_name, 1);
edge_n = size(edge_list, 1);

% check
for ei = 1 : edge_n
    i = edge_list(ei, 1);
    j = edge_list(ei, 2);
    if i == j, fprintf('Index pair error.\n'); pause; end
end

% load
I = cell(img_n, 1);
for i = 1 : img_n
    I{i} = imread(in_name{i});
end
% preprocessing
% T{i} is the coordinate transform of image i that change from the image
% coordinate (origin at the upper-left corner of the image, x towards right, y towards down) to
% the standard coordinate (origin at the center of the image, x towards
% right, y towards up)
max_size = 1000 * 1000;
imgw = zeros(img_n, 1);
imgh = zeros(img_n, 1);
resized_in_name = cell(img_n, 1);
T = cell(img_n, 1);
for i = 1 : img_n
    if numel(I{i}(:, :, 1)) > max_size % downsample
        I{i} = imresize(I{i}, sqrt(max_size / numel(I{i}(:, :, 1))));
    end
    resized_in_name{i} = sprintf('%s_resized%s', in_name{i}(1:end-4), in_name{i}(end-3:end));
    imwrite(I{i}, resized_in_name{i});
    imgw(i) = size(I{i}, 2);
    imgh(i) = size(I{i}, 1);
    T{i} = [1  0 -(imgw(i)+1)/2; ...
            0 -1  (imgh(i)+1)/2; ...
            0  0  1];
end
for i = 1 : img_n
    I{i} = im2double(I{i});
end

%% compute homography of image pairs
H = cell(img_n, img_n); % H{i,j} is the transformation from I{j} to I{i} (using standard coordinates)
for i = 1 : img_n
    H{i, i} = eye(3);
end
for ei = 1 : edge_n
    i = edge_list(ei, 1);
    j = edge_list(ei, 2);
    tmpH = sift_mosaic(I{i}, I{j}); % tmpH is the transform from I{i} to I{j} (using image coordinates)
    H{j, i} = T{j} * tmpH * inv(T{i}); % change of coord.
    H{i, j} = inv(H{j, i});
end
mdltpara = [];
mdlt_results = [];

%% for every image I, compute the homography between I and I_ref, and the homo. between I and I_target
% construct graph
G = sparse([edge_list(:, 1); edge_list(:, 2)], ...
           [edge_list(:, 2); edge_list(:, 1)], ...
           [ones(edge_n, 1); ones(edge_n, 1)], img_n, img_n); % construct a graph represented by a sparse mat
% compute homography between ith img and REFERENCE img
for i = 1 : img_n
    if i == ref, continue; end
    [~, path, ~] = graphshortestpath(G, i, ref); % path from ith img to ref img
    tmpH = eye(3);
    for ppi = 1 : numel(path)-1
        tmpH = H{path(ppi+1), path(ppi)} * tmpH;
    end
    H{ref, i} = tmpH;
    H{i, ref} = inv(H{ref, i});
end
% compute homography between ith img and TARGET img
for i = 1 : img_n
    if i == tar, continue; end
    [~, path, ~] = graphshortestpath(G, i, tar); % path from ith img to tar img
    tmpH = eye(3);
    for ppi = 1 : numel(path)-1
        tmpH = H{path(ppi+1), path(ppi)} * tmpH;
    end
    H{tar, i} = tmpH;
    H{i, tar} = inv(H{tar, i});
end
% normalize all homography s.t. h9 = 1
for hi = 1 : numel(H)
    if H{hi}
        H{hi} = H{hi} / H{hi}(3, 3);
    end
end

%% let H0 be the homography from reference to target, extract its parameters for computing our warp
H0 = H{tar, ref};
A = H0(1:2, 1:2);
t = H0(1:2, 3);
h31 = H0(3, 1);
h32 = H0(3, 2);
B = A - t*[h31 h32];
c = sqrt(h31^2+h32^2);
theta = atan2(-h32, -h31); % for compute ub1 and ub2

%% compute ub1 and ub2
% compute cu
cu = zeros(img_n, 1);
for i = 1 : img_n
    [tmpx, tmpy] = apply_transform(0, 0, H{ref, i}); % transform center (0, 0) of ith img to ref img
    [cu(i), ~] = apply_transform(tmpx, tmpy, [cos(theta), sin(theta), 0; ...
                                             -sin(theta), cos(theta), 0; ...
                                                       0,          0, 1]);
end
oriub1 = min(cu);
oriub2 = max(cu);
% ub1 = oriub1; ub2 = oriub2;

s_itv = 20; % sample interval
[offset_table1, offset_table2] = meshgrid(-300:10:600, -300:10:200);
totalcost_table = zeros(size(offset_table1));
for oi = 1 : size(offset_table1, 1)
for oj = 1 : size(offset_table1, 2)
    ub1 = oriub1 + offset_table1(oi, oj);
    ub2 = oriub2 + offset_table2(oi, oj);
    if ub2 - ub1 < 120
        totalcost_table(oi, oj) = nan;
        continue;
    end
    c1para = compute_c1_warp_coeff(H0, t, c, theta, ub1, ub2, zeroR_ON);
    
    % Jacobian of each image
    invR = [cos(theta), sin(theta), 0; ...
           -sin(theta), cos(theta), 0; ...
                     0,          0, 1];
    cost_list = zeros(1, img_n);
    for i = 1 : img_n
        x = linspace(1, imgw(i), ceil(imgw(i)/s_itv));
        y = linspace(1, imgh(i), ceil(imgh(i)/s_itv));
        [x, y] = meshgrid(x, y);
        [x, y] = apply_transform(x, y, T{i}); % coord change
        tmpH = invR * H{ref, i}; % tmpH maps (x,y) in ith img to (u,v) in ref img
        [u, v] = apply_transform(x, y, tmpH);
        reg1_mask = u < ub1;
        reg2_mask = u > ub1 & u < ub2;
        reg3_mask = u > ub2;

        J11_map = zeros(size(x));
        J12_map = zeros(size(x));
        J21_map = zeros(size(x));
        J22_map = zeros(size(x));
        for reg = 1 : 3 % for each region
            if reg == 1,        reg_mask = u < ub1;
            elseif reg == 3,    reg_mask = u > ub2;
            else                reg_mask = ~(u < ub1 | u > ub2);
            end
            if nnz(reg_mask) == 0, continue; end
            x1 = x(reg_mask);
            y1 = y(reg_mask);
            u1 = u(reg_mask);
            v1 = v(reg_mask);

            if reg == 1 % region 1
                A = c1para.H * H{ref, i};
            elseif reg == 2 % region 2
                A = tmpH; % as described above
            elseif reg == 3 % region 3
                A = c1para.S * H{ref, i};
            end
            A = A / A(3, 3);
            h1 = A(1, 1); h2 = A(1, 2); h3 = A(1, 3); h4 = A(2, 1); h5 = A(2, 2); h6 = A(2, 3); h7 = A(3, 1); h8 = A(3, 2);

            J11 = h1./(h7*x1 + h8*y1 + 1) - (h7*(h3 + h1*x1 + h2*y1))./(h7*x1 + h8*y1 + 1).^2;  J12 = h2./(h7*x1 + h8*y1 + 1) - (h8*(h3 + h1*x1 + h2*y1))./(h7*x1 + h8*y1 + 1).^2;
            J21 = h4./(h7*x1 + h8*y1 + 1) - (h7*(h6 + h4*x1 + h5*y1))./(h7*x1 + h8*y1 + 1).^2;  J22 = h5./(h7*x1 + h8*y1 + 1) - (h8*(h6 + h4*x1 + h5*y1))./(h7*x1 + h8*y1 + 1).^2;
            if reg == 2
                if c1para.zeroR_ON == 0
                    JT11 = (2*c1para.a(1)*u1 + c1para.a(2)) .* v1 + (2*c1para.b(1)*u1 + c1para.b(2));
                    JT12 = c1para.a(1)*u1.^2 + c1para.a(2)*u1 + c1para.a(3);
                    JT21 = (2*c1para.e(1)*u1 + c1para.e(2)) .* v1 + (2*c1para.f(1)*u1 + c1para.f(2));
                    JT22 = c1para.e(1)*u1.^2 + c1para.e(2)*u1 + c1para.e(3);
                else
                    JT11 = (3*c1para.a(1)*u1.^2 + 2*c1para.a(2)*u1 + c1para.a(3)) .* v1 + (2*c1para.b(1)*u1 + c1para.b(2));
                    JT12 = c1para.a(1)*u1.^3 + c1para.a(2)*u1.^2 + c1para.a(3)*u1 + c1para.a(4);
                    JT21 = (3*c1para.e(1)*u1.^2 + 2*c1para.e(2)*u1 + c1para.e(3)) .* v1 + (2*c1para.f(1)*u1 + c1para.f(2));
                    JT22 = c1para.e(1)*u1.^3 + c1para.e(2)*u1.^2 + c1para.e(3)*u1 + c1para.e(4);
                end
                tmp11 = JT11.*J11 + JT12.*J21;
                tmp12 = JT11.*J12 + JT12.*J22;
                tmp21 = JT21.*J11 + JT22.*J21;
                tmp22 = JT21.*J12 + JT22.*J22;
                J11 = tmp11; J12 = tmp12; J21 = tmp21; J22 = tmp22;
            end
            J11_map(reg_mask) = J11;
            J12_map(reg_mask) = J12;
            J21_map(reg_mask) = J21;
            J22_map(reg_mask) = J22;
        end
        % As-GlobalSimilar-As-Possible cost
        avg_alpha = (sum(J11_map(:)) + sum(J22_map(:))) / (numel(x)*2);
        avg_beta = (sum(J21_map(:)) + (-sum(J12_map(:)))) / (numel(x)*2);
        cost = sum( (J11_map(:)-avg_alpha).^2 + (J12_map(:)-(-avg_beta)).^2 + ...
                    (J21_map(:)-avg_beta).^2  + (J22_map(:)-avg_alpha).^2 );
        % done
        cost_list(i) = cost;
    end
    totalcost = sum(cost_list);
    totalcost_table(oi, oj) = totalcost;
end
end
% figure; surf(totalcost_table);
% fprintf('min cost = %f\n', min(totalcost_table(:)));
[idx1,idx2] = find(totalcost_table==min(totalcost_table(:)));
ub1 = oriub1 + offset_table1(idx1, idx2);
ub2 = oriub2 + offset_table2(idx1, idx2);
% fprintf('init (u1,u2) = (%.0f,%.0f). offset = (%.0f, %.0f).\n', oriub1, oriub2, offset_table1(idx1, idx2), offset_table2(idx1, idx2));

%% Compute parameters/coefficients of our warp
c1para = compute_c1_warp_coeff(H0, t, c, theta, ub1, ub2, zeroR_ON);

%% c1 warp for each image
img_grid_size = 10; % for warping
[c1out, c1omask, Td] = texture_mapping(resized_in_name, imgw, imgh, img_grid_size, T, warp_type, H, ref, tar, c1para, mdltpara, 0, 'white');
for i = 1 : img_n
    %figure; imshow(c1out{i}); % the warped result of each image
end

for i = 1 : img_n
    eval(['delete ' resized_in_name{i}]);
end

%% composite by linear blending
for i = 1 : img_n
    if i == 1
        out = c1out{1};
        out_mask = c1omask{1};
        % center of out
        [r, c] = find(c1omask{1}(:, :, 1));
        out_center = [mean(r) mean(c)];
    else % blend out and c1out{i}
        % center of c1out{i}
        [r, c] = find(c1omask{i}(:, :, 1));
        out_i_center = [mean(r) mean(c)];
        % compute weighted mask
        vec = out_i_center - out_center; % vector from out_center to out_i_center
        intsct_mask = c1omask{i}(:, :, 1) & out_mask(:, :, 1); % 1 channel
        out_only_mask = out_mask(:, :, 1) - intsct_mask;
        out_i_only_mask = c1omask{i}(:, :, 1) - intsct_mask;
        
        [r, c] = find(intsct_mask(:, :, 1));
        idx = sub2ind(size(c1omask{i}(:, :, 1)), r, c);
        out_wmask = zeros(size(c1omask{i}(:, :, 1)));
        proj_val = (r - out_center(1))*vec(1) + (c- out_center(2))*vec(2); % inner product
        out_wmask(idx) = (proj_val - (min(proj_val)+(1e-3))) / ...
                         ((max(proj_val)-(1e-3)) - (min(proj_val)+(1e-3))); % weight map (of overlapped area) for c1out{i}, 1 channel
        % blending
        mask1 = out_mask(:, :, 1)&(out_wmask==0);
        mask2 = out_wmask;
        mask3 = c1omask{i}(:, :, 1)&(out_wmask==0);
        mask1 = cat(3, mask1, mask1, mask1); mask2 = cat(3, mask2, mask2, mask2); mask3 = cat(3, mask3, mask3, mask3);
        out = out.*(mask1+(1-mask2).*(mask2~=0)) + c1out{i}.*(mask2+mask3);
        % update
        out_mask = out_mask | c1omask{i};
        out_center = out_i_center; % update out_center by assign center of c1out{i}
    end
end
c1all = out;
denom = zeros(size(c1out{1}));
for i = 1 : img_n
    denom = denom + c1omask{i};
end
bgcolor = 'white'; % background color
if strcmp(bgcolor, 'white')
    c1all(denom==0) = 1;
else
    c1all(denom==0) = 0;
end
figure; imshow(c1all); % the final result
rename = sprintf('%sresult%s', in_name{1}(1:end-4), in_name{1}(end-3:end));
imwrite(c1all, rename);