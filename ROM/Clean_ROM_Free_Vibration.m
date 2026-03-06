%% Reduced Order Model (ROM)
%% 2 DOF System
clear
dof = 2;

%% Beam Inputs
L   = 215e-3;      % m
t   = 0.5e-3;      % m
b   = 7.5e-3;      % m
E   = 4e9;         % Pa
rho = 1240;        % kg m^-3

A = b * t;
I = b * t^3 / 12;

%% Arch Geometry
beta_beam   = 0;
c_beam      = 0.8;
num_nodes   = 100;
initial_disp = 0;

[x,y_0,~] = discretise_arch(beta_beam,c_beam,L,num_nodes);
y_0 = y_0';

%% Mode Shapes

phi_b      = zeros(length(x),dof);
phi_b_norm = zeros(length(x),dof);

lambda_i = zeros(1,dof);
delta_i  = zeros(1,dof);

ig = 1;

trnsc_eq = @(l) cos(l).*cosh(l) - 1;

for n = 1:dof

    lambda_i(n) = fzero(trnsc_eq,ig);
    ig = lambda_i(n) + 2;

    delta_i(n) = (cosh(lambda_i(n)) - cos(lambda_i(n))) ...
               /(sinh(lambda_i(n)) - sin(lambda_i(n)));

end

for n = 1:dof

    frac = (lambda_i(n) * x) / L;

    phi_b(:,n) = cosh(frac) - cos(frac) ...
               - delta_i(n) * (sinh(frac) - sin(frac));

    phi_b_norm(:,n) = phi_b(:,n) / max(phi_b(:,n));

end

% Modal coefficients for initial arch shape
gamma = phi_b_norm \ y_0;

%% Two–Mode Projection of Arch Shape

y_recon = phi_b_norm * gamma;
residual = y_0 - y_recon;

figure
plot(x,y_0,'k','DisplayName',"Original")
hold on
plot(x,y_recon,'b','DisplayName',"Reconstructed (Two Mode)")
plot(x,residual,'--r','DisplayName',"Residual")

grid on
legend

xlabel("x")
ylabel("y")

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

        M(i,j) = rho * A * trapz(x,integrand_M);
        G(i,j) = trapz(x,integrand_G);
        K(i,j) = E * I * trapz(x,integrand_K);

    end
end

%% Equation of Motion Parameters

P = 0; % Synthesised-in-shape arch - No pre-stress

M11 = M(1,1);
M22 = M(2,2);

K11 = K(1,1);
K22 = K(2,2);

G11 = G(1,1);
G22 = G(2,2);

%% Equilibrium Solution

dVda = @(a) [

    K11*(a(1)-gamma(1)) ...
    - G11*a(1)*(P + ((A*E)/(2*L))*(G11*gamma(1)^2 + G22*gamma(2)^2)) ...
    + ((A*E)/(2*L))*G11*a(1)*(G11*a(1)^2 + G22*a(2)^2);

    K22*(a(2)-gamma(2)) ...
    - G22*a(2)*(P + ((A*E)/(2*L))*(G11*gamma(1)^2 + G22*gamma(2)^2)) ...
    + ((A*E)/(2*L))*G22*a(2)*(G11*a(1)^2 + G22*a(2)^2)

];

igu1 = [-1e3 ; -1e3];

a_1 = fsolve(dVda,igu1);

%% Equilibrium Shape

equilibria_Y = a_1(1)*phi_b_norm(:,1) ...
             + a_1(2)*phi_b_norm(:,2);

%% Initial State Vector

a1_0 = a_1(1) + initial_disp;
a2_0 = a_1(2);

a1_dot_0 = 0;
a2_dot_0 = 0;

a_vec_0 = [a1_0 ; a1_dot_0 ; a2_0 ; a2_dot_0];

%% Time Span

tspan = [0 2];

%% Solve EOM

[t,a_vec] = ode45(@(t,a_vec) ...
    archEOM(t,a_vec,M11,M22,K11,K22,G11,G22,gamma(1),gamma(2),P,A,E,L), ...
    tspan,a_vec_0);

%% Reconstruct Displacement Field

y = phi_b_norm(:,1)*a_vec(:,1)' ...
  + phi_b_norm(:,2)*a_vec(:,3)';

%% Time History Plot

figure

POI = round(num_nodes/2);

plot(t,y(POI,:))
grid on

xlabel("Time / s")
ylabel("Displacement / m")

%% Spatiotemporal Surface

figure

N = size(t);
step = ceil(N/100);
idx = 1:step:N;

y_surf = y(:,idx);
t_surf = t(idx,:);

surf(x,t_surf,y_surf')

shading interp
colormap turbo
axis tight

xlabel('Beam Length (m)')
ylabel('Time (s)')
zlabel('Deflection (m)')
title('Spatiotemporal Surface')

view(45,45)

%% Potential Energy Surface

a1_range = linspace(-5e-3,5e-3,200);
a2_range = linspace(-5e-3,5e-3,200);

[A1,A2] = meshgrid(a1_range,a2_range);

V = 0.5*K11.*(A1-2*gamma(1)).*A1 ...
  + 0.5*K22.*(A2-2*gamma(2)).*A2 ...
  - 0.5*(P + (A*E/(2*L))*(G11*gamma(1)^2 + G22*gamma(2)^2)) ...
      .* (G11*A1.^2 + G22*A2.^2) ...
  + (A*E/(8*L))*(G11*A1.^2 + G22*A2.^2).^2;

figure
surf(A1,A2,V)

shading interp
colormap turbo
view(45,35)

xlabel('a_1')
ylabel('a_2')
zlabel('Potential Energy')

title('Potential Energy Surface')

%% ------------------------------------------------------------------------
%% FUNCTIONS
%% ------------------------------------------------------------------------

function dydt = archEOM(~,a_vec,M11,M22,K11,K22,G11,G22,gamma1,gamma2,P,A,E,L)

a1     = a_vec(1);
a1_dot = a_vec(2);

a2     = a_vec(3);
a2_dot = a_vec(4);

% Non-linear terms
axial_term = P + (A*E/(2*L))*(G11*gamma1^2 + G22*gamma2^2);
cubic_term = (A*E/(2*L))*(G11*a1^2 + G22*a2^2);

% Accelerations
a1_dd = -(1/M11)*( ...
      K11*(a1-gamma1) ...
    - G11*a1*axial_term ...
    + G11*a1*cubic_term );

a2_dd = -(1/M22)*( ...
      K22*(a2-gamma2) ...
    - G22*a2*axial_term ...
    + G22*a2*cubic_term );

% First order system
dydt = zeros(4,1);

dydt(1) = a1_dot;
dydt(2) = a1_dd;

dydt(3) = a2_dot;
dydt(4) = a2_dd;

end

function [X,Y,theta] = discretise_arch(Beta_A,c_A,L_beam,num_points)

% Models an arch using the formulation from L.N. Virgin
% Geometry is then discretised using the specified number of points

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
