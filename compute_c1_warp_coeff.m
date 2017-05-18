function c1para = compute_c1_warp_coeff(H, t, c, theta, ub1, ub2, zeroR_ON)
% c1para has the following attributes:
% theta, H, ub1, ub2, a, b, e, f, S

if ~exist('zeroR_ON', 'var'), zeroR_ON = 0; end

%% homography matrix
cos_theta = cos(theta);
sin_theta = sin(theta);
%
c1para.theta = theta;
c1para.H = H;
c1para.zeroR_ON = zeroR_ON;

if ~exist('ub1', 'var') || ~exist('ub1', 'var')    
    c1para.ub1 = inf;
    c1para.ub2 = inf;
    c1para.a = NaN;
    c1para.b = NaN;
    c1para.e = NaN;
    c1para.f = NaN;
    c1para.S = NaN(3,3);
    return
end

p1 = [H(1, 1) H(1, 2)] * [cos_theta; sin_theta];
p2 = [H(2, 1) H(2, 2)] * [cos_theta; sin_theta];
q1 = [H(1, 1) H(1, 2)] * [-sin_theta; cos_theta];
q2 = [H(2, 1) H(2, 2)] * [-sin_theta; cos_theta];

%% interpolation coeff
%% A and E
% A(u)
% para. of A(u)
G = [ub1^2, ub1, 1, 0; ...
     ub2^2, ub2, 1, 1; ...
     2*ub1, 1, 0, 0; ...
     2*ub2, 1, 0, 0];
g = [q1/(1-c*ub1); 0; c*q1/(1-c*ub1)^2; 0];
sol = inv(G)*g;
a = sol(1:3);
beta = sol(4);
% E(u)
% para. of e(u)
G = [ub1^2, ub1, 1, 0; ...
     ub2^2, ub2, 1, 1; ...
     2*ub1, 1, 0, 0; ...
     2*ub2, 1, 0, 0];
g = [q2/(1-c*ub1); 0; c*q2/(1-c*ub1)^2; 0];
sol = inv(G)*g;
e = sol(1:3);
alpha = -sol(4);
%% compute again if angle of R(x,y) is constrained to be zero
scaling = sqrt(alpha^2 + beta^2);
if zeroR_ON
    % A(u)
    G = [ub1^3, ub1^2, ub1, 1; ...
         ub2^3, ub2^2, ub2, 1; ...
         3*ub1^2, 2*ub1, 1, 0; ...
         3*ub2^2, 2*ub2, 1, 0];
    g = [q1/(1-c*ub1); -scaling*sin(theta); c*q1/(1-c*ub1)^2; 0];
    sol = inv(G)*g;
    a = sol;
    beta = scaling*sin(theta);
    % E(u)
    G = [ub1^3, ub1^2, ub1, 1; ...
         ub2^3, ub2^2, ub2, 1; ...
         3*ub1^2, 2*ub1, 1, 0; ...
         3*ub2^2, 2*ub2, 1, 0];
    g = [q2/(1-c*ub1); scaling*cos(theta); c*q2/(1-c*ub1)^2; 0];
    sol = inv(G)*g;
    e = sol;
    alpha = scaling*cos(theta);
end
%% B and F
% B(u)
% para. of B(u) and similarity transform
G = [ub1^2, ub1, 1, 0; ...
     ub2^2, ub2, 1, 1; ...
     2*ub1, 1, 0, 0; ...
     2*ub2, 1, 0, 0];
g = [(p1*ub1+H(1, 3))/(1-c*ub1); ...
     alpha*ub2; ...
     (p1+c*H(1, 3))/(1-c*ub1)^2; ...
     alpha];
sol = inv(G)*g;
b = sol(1:3);
transx = -sol(4);
% F(u)
% para. of F(u) and similarity transform
G = [ub1^2, ub1, 1, 0; ...
     ub2^2, ub2, 1, 1; ...
     2*ub1, 1, 0, 0; ...
     2*ub2, 1, 0, 0]; % the same w/ G of B(u)
g = [(p2*ub1+H(2, 3))/(1-c*ub1); ...
     beta*ub2; ...
     (p2+c*H(2, 3))/(1-c*ub1)^2; ...
     beta];
sol = inv(G)*g;
f = sol(1:3);
transy = -sol(4);

%% similarity transform in terms of (x,y)
S = [alpha, -beta, transx; ...
      beta, alpha, transy; ...
         0,     0,      1]; % this is in terms of (u,v)
S = S * [ cos_theta, sin_theta, 0; ...
         -sin_theta, cos_theta, 0; ...
                  0,         0, 1];
if zeroR_ON
    S(2,1) = 0;
    S(1,2) = 0; % clean numerical nonzero
end

c1para.ub1 = ub1;
c1para.ub2 = ub2;
c1para.a = a;
c1para.b = b;
c1para.e = e;
c1para.f = f;
c1para.S = S;