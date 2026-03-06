%% Reduced Order Model – Linearised Natural Frequencies
%% 2 DOF System

clear
dof = 2;

%% Beam Inputs

L   = 150e-3;     % m
t   = 0.8e-3;     % m
b   = 7.5e-3;     % m
E   = 4e9;        % Pa
rho = 1240;       % kg m^-3

A = b * t;
I = b * t^3 / 12;

%% Arch Geometry

beta_beam = 0;
c_beam    = 0.8;
L_beam    = L;

num_nodes = 100;

[x,y,~] = discretise_arch(beta_beam,c_beam,L_beam,num_nodes);

x = x';
y = y';

%% Mode Shapes

phi_b      = zeros(length(x),dof);
phi_b_norm = zeros(length(x),dof);

lambda_i = zeros(1,dof);
delta_i  = zeros(1,dof);

ig = 1;

trnsc_eq = @(l) cos(l)*cosh(l) - 1;

for n = 1:dof

    lambda_i(n) = fzero(trnsc_eq,ig);
    ig = lambda_i(n) + 2;

    delta_i(n) = (cosh(lambda_i(n)) - cos(lambda_i(n))) ...
               /(sinh(lambda_i(n)) - sin(lambda_i(n)));

end

for n = 1:dof

    frac = (lambda_i(n)*x)/L;

    phi_b(:,n) = cosh(frac) - cos(frac) ...
               - delta_i(n)*(sinh(frac) - sin(frac));

    phi_b_norm(:,n) = phi_b(:,n)/max(phi_b(:,n));

end

%% Matrix Formulation

dx = x(2) - x(1);

dphi  = zeros(size(phi_b_norm));
ddphi = zeros(size(phi_b_norm));

for n = 1:dof
    dphi(:,n)  = gradient(phi_b_norm(:,n),dx);
    ddphi(:,n) = gradient(dphi(:,n),dx);
end

M = zeros(dof,dof);
G = zeros(dof,dof);
K = zeros(dof,dof);

for i = 1:dof
    for j = 1:dof

        integrand_M = phi_b_norm(:,i) .* phi_b_norm(:,j);
        integrand_G = dphi(:,i)       .* dphi(:,j);
        integrand_K = ddphi(:,i)      .* ddphi(:,j);

        M(i,j) = rho*A*trapz(x,integrand_M);
        G(i,j) = trapz(x,integrand_G);
        K(i,j) = E*I*trapz(x,integrand_K);

    end
end

%% Equilibrium Parameters

gamma = zeros(dof,1);

gamma(1) = 0.082e-3;
gamma(2) = -0.077e-3;

P = 184.1;

K1 = K(1,1);
K2 = K(2,2);

G1 = G(1,1);
G2 = G(2,2);

%% Solve for Equilibrium Coordinates

dVda = @(a) [

    K1*(a(1)-gamma(1)) ...
    - G1*a(1)*(P + ((A*E)/(2*L))*(G1*gamma(1)^2 + G2*gamma(2)^2)) ...
    + ((A*E)/(2*L))*G1*a(1)*(G1*a(1)^2 + G2*a(2)^2);

    K2*(a(2)-gamma(2)) ...
    - G2*a(2)*(P + ((A*E)/(2*L))*(G1*gamma(1)^2 + G2*gamma(2)^2)) ...
    + ((A*E)/(2*L))*G2*a(2)*(G1*a(1)^2 + G2*a(2)^2)

];

igu1 = [-1e3 ; -1e3];

a_1 = fsolve(dVda,igu1);

%% Linearised Natural Frequencies

a  = a_1;
gt = gamma';
at = a';

Gt = G';

Kbar = K ...
     - (P + ((A*E)/(2*L))*(gt*G*gamma - at*G*a))*G ...
     + ((A*E)/L)*G*a*at*Gt;

Amat = inv(M)*Kbar;

lambda = eig(Amat);

wn = sqrt(lambda);
fn = wn/(2*pi);

disp(fn)

%% Stability Check

chk = isreal(wn);

if chk == [1,1]

    if y((length(x)/2)) > 0

        disp("Up")
        disp("Stable")

        diff1 = abs(115 - fn(1));
        diff2 = abs(173 - fn(2));

        per_diff1 = (diff1/115)*100;
        per_diff2 = (diff2/173)*100;

        fprintf("Percentage Difference %d for Mode 1 and %d for Mode 2", ...
            per_diff1,per_diff2)

    else

        disp("Down")
        disp("Stable")

        diff1 = abs(102 - fn(1));
        diff2 = abs(171 - fn(2));

        per_diff1 = (diff1/102)*100;
        per_diff2 = (diff2/171)*100;

        fprintf("Percentage Difference %d for Mode 1 and %d for Mode 2", ...
            per_diff1,per_diff2)

    end

else

    disp("Unstable")

end

%% Arch Geometry Plot

figure

x_mm = 1000*x;
y_mm = 1000*y;

plot(x_mm,y_mm,linewidth=2)

xlabel("x (mm)")
ylabel("y (mm)")

grid on

%% ------------------------------------------------------------------------
%% FUNCTIONS
%% ------------------------------------------------------------------------

function [X,Y,theta] = discretise_arch(Beta_A,c_A,L_beam,num_points)

% Models an arch using the formulation from L.N. Virgin
% Geometry is discretised using the specified number of points

N = num_points;

s = linspace(0,1,N);

theta = (1-2*s).*Beta_A ...
      + c_A.*s.*(1-s).*(1-2*s);

X_nd = cumtrapz(s,cos(theta));
Y_nd = cumtrapz(s,sin(theta));

scale = L_beam / X_nd(end);

X = X_nd * scale;
Y = Y_nd * scale;

end