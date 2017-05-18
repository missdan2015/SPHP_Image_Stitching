function [out_x, out_y] = c1_warp_ver2(x, y, c1para)

theta = c1para.theta;
H = c1para.H;
if c1para.ub1 == inf && c1para.ub2 == inf % pure homography
    out_x = (H(1, 1)*x + H(1, 2)*y + H(1, 3)) ./ ...
            (H(3, 1)*x + H(3, 2)*y + H(3, 3));
    out_y = (H(2, 1)*x + H(2, 2)*y + H(2, 3)) ./ ...
            (H(3, 1)*x + H(3, 2)*y + H(3, 3));
    return
end
a = c1para.a;
b = c1para.b;
e = c1para.e;
f = c1para.f;

u =  cos(theta)*x + sin(theta)*y;
v = -sin(theta)*x + cos(theta)*y;
% mask
H_mask = u <= c1para.ub1; % homography
Interp_mask = u <= c1para.ub2 & u > c1para.ub1;
S_mask = u > c1para.ub2; % similarity

% homography transform
x1 = x(H_mask);
y1 = y(H_mask);
outx1 = (H(1, 1)*x1 + H(1, 2)*y1 + H(1, 3)) ./ ...
        (H(3, 1)*x1 + H(3, 2)*y1 + H(3, 3));
outy1 = (H(2, 1)*x1 + H(2, 2)*y1 + H(2, 3)) ./ ...
        (H(3, 1)*x1 + H(3, 2)*y1 + H(3, 3));
% similarity transform
x3 = x(S_mask);
y3 = y(S_mask);
S = c1para.S;
outx3 = S(1, 1)*x3 + S(1, 2)*y3 + S(1, 3);
outy3 = S(2, 1)*x3 + S(2, 2)*y3 + S(2, 3);
% interpolation
u2 = u(Interp_mask);
v2 = v(Interp_mask);
if c1para.zeroR_ON == 0
    outx2 = (a(1)*u2.^2 + a(2)*u2 + a(3)) .* v2 + (b(1)*u2.^2 + b(2)*u2 + b(3));
    outy2 = (e(1)*u2.^2 + e(2)*u2 + e(3)) .* v2 + (f(1)*u2.^2 + f(2)*u2 + f(3));
else
    outx2 = (a(1)*u2.^3 + a(2)*u2.^2 + a(3)*u2 + a(4)) .* v2 + (b(1)*u2.^2 + b(2)*u2 + b(3));
    outy2 = (e(1)*u2.^3 + e(2)*u2.^2 + e(3)*u2 + e(4)) .* v2 + (f(1)*u2.^2 + f(2)*u2 + f(3));
end

out_x = zeros(size(x));
out_x(H_mask) = outx1;
out_x(Interp_mask) = outx2;
out_x(S_mask) = outx3;

out_y = zeros(size(x));
out_y(H_mask) = outy1;
out_y(Interp_mask) = outy2;
out_y(S_mask) = outy3;