%% Corotational FEM – Forced Vibration
%% Arched Beam

fprintf("Everything saved?")
pause

clear all; clc;
% close all;

%% Beam Inputs

L_beam = 215e-3;                    % Length (m)
b = 7.5e-3;                         % Width (m)
t = 0.5e-3;                         % Thickness (m)
E = 2636e6;                         % Young's modulus (Pa)
rho = 1205;                         % Density (kg m^-3)

I = b * t^3 / 12;
A = b * t;

%% Arch Geometry Inputs

c_beam = 0.8;                       % L.N. Virgin equation
beta_beam = 0;                      % L.N. Virgin equation

%% Simulation Inputs

forced_vibration = "Mono-harmonic"; % "Mono-harmonic" or "Swept-sine"

num_e = 24;
sim_time = 30;                      % Simulation time (s)
damp_fac = 0.8;                     % Damping factor
accel_magnitude = 5.295;            % Base acceleration amplitude (g)
force_freq = 39;                    % Mono-harmonic forcing frequency (Hz)
time_scale = 100;                   % Percentage of critical timestep

f_low = 10;                         % Starting sweep frequency (Hz)
R = 4;                              % Sweep rate (oct/min)

measurement_point = 15;             % cm from built-in end

%% Discretisation

nodal_dof = 3;
num_nodes = num_e + 1;
num_dof = nodal_dof * num_nodes;
mid_point = (num_e / 2) + 1;

%% Geometry

[x_o, y_o, th_o] = discretise_arch(beta_beam,c_beam,L_beam,num_nodes);
x_o = x_o';
y_o = y_o';
th_o = th_o';

%% Initial Element Lengths

l_o = zeros(num_e,1);

for n = 1:num_e
    dx = x_o(n+1) - x_o(n);
    dy = y_o(n+1) - y_o(n);
    l_o(n) = sqrt(dx^2 + dy^2);
end

%% Simulation Controls

c_wave = sqrt(E/rho);
crit_timestep = min(l_o)/c_wave;
fprintf('Critical timestep %d \n',crit_timestep)

timestep = (time_scale/100) * crit_timestep;
num_itr = ceil(sim_time / timestep);
fprintf('% d iterations required, press any key to progress or ctrl + c to cancel \n',num_itr)
pause

if timestep > crit_timestep
    fprintf('Timestep too large (%d > %d)',timestep,crit_timestep)
    return
end

%% System Construction

% Lumped mass matrix as described in Wiebe PhD.
[Me] = corotational_local_mass_matrix(num_e,nodal_dof,A,l_o,rho);
[Mm] = formulate_corotational_global_mass_matrix(num_e,num_dof,Me);

invM = inv(Mm);                     % Pre-calculate inverse mass

%% Initialisation

d_o     = zeros(num_dof,1);         % Initial displacement vector
d_dot_o = zeros(num_dof,1);         % Initial velocity vector

% Pre-allocate storage.
y_itr     = zeros(num_nodes,num_itr+1);
y_dot_itr = zeros(num_nodes,num_itr+1);

% Set initial history.
x_itr_o = x_o + d_o(1:3:end);
y_itr(:,1) = y_o + d_o(2:3:end);
th_itr_o = th_o + d_o(3:3:end);

y_dot_itr(:,1) = 0;

% State vectors.
u_curr = d_o;
v_curr = d_dot_o;

%% Co-rotational Initialisation

[Beta,Le] = calculate_length_n_beta(num_e,x_itr_o,y_itr(:,1));
[B] = corotational_B_matrix(num_e,Beta,Le);

[N,M1,M2] = element_internal_forcing(num_e,E,A,I,x_o,y_o,th_o, ...
    x_itr_o,y_itr(:,1),th_itr_o);
[F_int] = global_internal_forcing(num_e,nodal_dof,B,N,M1,M2);

[Ke] = corotational_local_stiffness_matrix(num_e,nodal_dof,E,A,I, ...
    l_o,Le,B,Beta,N,M1,M2);
[Kk] = formulate_global_stiffness_matrix(num_e,num_dof,Ke);

%% Time Stepping Inputs

fprintf('Starting simulation...\n');
time_vec = (0:num_itr)'*timestep;

%% Base Excitation

if forced_vibration == "Mono-harmonic"

    base_accel = accel_magnitude * 9.81 .* sin(2*pi*force_freq.*time_vec);

elseif forced_vibration == "Swept-sine"

    alpha = (R/60) * log(2);
    phi = (2*pi*f_low/alpha) * (exp(alpha*time_vec)-1);

    base_accel = accel_magnitude * 9.81 .* sin(phi);

end

r_vec = zeros(num_dof,1);
r_vec(2:3:end) = 1;

%% Time Stepping Loop: Explicit Euler

for itr = 1:num_itr

    % Update geometry based on current displacement.
    x_curr = x_o + u_curr(1:3:end);
    y_curr = y_o + u_curr(2:3:end);
    th_curr = th_o + u_curr(3:3:end);

    if max(abs(u_curr(2:3:end))) > 100*max(y_o)
        fprintf('Warning: divergence in iteration %d',itr)
        return
    end

    % Index relevant force vector.
    F_ext = -Mm * r_vec * base_accel(itr);

    % Mass-proportional damping.
    F_damp = damp_fac * Mm * v_curr;

    % Equation of motion.
    accel = invM * (F_ext - F_int - F_damp);

    accel(1:3) = 0;
    u_curr(1:3) = 0;
    v_curr(1:3) = 0;

    accel(end-2:end) = 0;
    u_curr(end-2:end) = 0;
    v_curr(end-2:end) = 0;

    % Explicit time update.
    v_next = v_curr + timestep * accel;
    u_next = u_curr + timestep * v_next;

    % Store and advance.
    x_itr_1 = x_o + u_next(1:3:end);
    y_itr(:,itr+1) = y_o + u_next(2:3:end);
    th_itr_1 = th_o + u_next(3:3:end);
    y_dot_itr(:,itr+1) = v_next(2:3:end);

    u_curr = u_next;
    v_curr = v_next;

    [Beta,Le] = calculate_length_n_beta(num_e,x_itr_1,y_itr(:,itr+1));
    [B] = corotational_B_matrix(num_e,Beta,Le);
    [N,M1,M2] = element_internal_forcing(num_e,E,A,I,x_o,y_o,th_o, ...
        x_itr_1,y_itr(:,itr+1),th_itr_1);
    [F_int] = global_internal_forcing(num_e,nodal_dof,B,N,M1,M2);

    if mod(itr,round(num_itr/20)) == 0
        fprintf('Progress: %.0f%%\n',(itr/num_itr)*100);
    end

end

%% Time History Plot

figure

diff_x = x_o - measurement_point*1e-2;
[~,x_idx] = min(abs(diff_x));

subplot(2,1,1)
plot(time_vec,y_itr(x_idx,:))
xlabel('Time (s)')
ylabel('Deflection (m)')

subplot(2,1,2)
plot(time_vec,y_dot_itr(x_idx,:))
xlabel('Time (s)')
ylabel('Velocity (m/s)')

sgtitle("Time Series")

%% Outputs

y_numerical = y_itr(x_idx,:);
y_dot_numerical = y_dot_itr(x_idx,:);
base_accel_numerical = base_accel;
t_numerical = time_vec;

%% ------------------------------------------------------------------------
%% Functions
%% ------------------------------------------------------------------------

%% Bernoulli-Euler Beam: Lumped Mass Matrix, E.6 Wiebe PhD

function [Me] = corotational_local_mass_matrix(num_e,nodal_dof,A,L,rho)

    e_dof = nodal_dof * 2;
    Me = zeros(e_dof,e_dof,num_e);

    for n = 1:num_e
        Lw = L(n);
        const = 0.5 * (rho * A * Lw);
        Lw_12 = (Lw^2)/12;
        mat = diag([1 1 Lw_12 1 1 Lw_12]);
        Me(:,:,n) = const * mat;
    end

end

function [Mm] = formulate_corotational_global_mass_matrix(num_e,num_dof,Me)

    Mm = zeros(num_dof,num_dof);

    for n = 1:num_e

        % No rotation is applied to the lumped mass matrix.
        Me_w = Me(:,:,n);

        node1 = n;
        node2 = n+1;
        dof_map = [3*node1-2, 3*node1-1, 3*node1, ...
                   3*node2-2, 3*node2-1, 3*node2];

        Mm(dof_map,dof_map) = Mm(dof_map,dof_map) + Me_w;

    end

end

function [Beta,Le] = calculate_length_n_beta(num_e,x,y)

    Beta = zeros(num_e,1);
    Le = zeros(num_e,1);

    for n = 1:num_e
        dX = x(n+1) - x(n);
        dY = y(n+1) - y(n);

        Le(n) = sqrt(dX^2 + dY^2);
        Beta(n) = atan2(dY,dX);
    end

end

function [N,M1,M2] = element_internal_forcing(num_e,E,A,I, ...
    x_o,y_o,th_o,x_new,y_new,th_new)

    N = zeros(num_e,1);
    M1 = zeros(num_e,1);
    M2 = zeros(num_e,1);

    for n = 1:num_e

        dx_o = x_o(n+1) - x_o(n);
        dy_o = y_o(n+1) - y_o(n);
        beta_o = atan2(dy_o,dx_o);

        dx = x_new(n+1) - x_new(n);
        dy = y_new(n+1) - y_new(n);
        beta = atan2(dy,dx);

        Le_o = sqrt(dx_o^2 + dy_o^2);
        Le = sqrt(dx^2 + dy^2);

        c_o = cos(beta_o);
        s_o = sin(beta_o);
        c = cos(beta);
        s = sin(beta);

        sin_alpha = c_o*s - s_o*c;
        cos_alpha = c_o*c + s_o*s;
        alpha = atan2(sin_alpha,cos_alpha);

        theta_1c = (th_new(n) - th_o(n)) - alpha;
        theta_2c = (th_new(n+1) - th_o(n+1)) - alpha;

        EA_Lo = (E*A)/Le_o;
        EI_Lo = (E*I)/Le_o;

        N(n) = EA_Lo * (Le - Le_o);
        M1(n) = 2 * EI_Lo * (2*theta_1c + theta_2c);
        M2(n) = 2 * EI_Lo * (theta_1c + 2*theta_2c);

    end

end

function [B] = corotational_B_matrix(num_e,Beta,Le)

    B = zeros(3,6,num_e);

    for n = 1:num_e

        Le_w = Le(n);
        c = cos(Beta(n));
        s = sin(Beta(n));

        mat2 = [-c*Le_w, -s*Le_w,   0, c*Le_w, s*Le_w, 0; ...
                -s,        c,     Le_w, s,       -c,    0; ...
                -s,        c,       0,  s,       -c,  Le_w];

        B(:,:,n) = (1/Le(n)) * mat2;

    end

end

%% Global Stiffness Matrix: Euler Beam Elements

function [Ke] = corotational_local_stiffness_matrix(num_e,nodal_dof,E,A,I, ...
    Le_o,Le,B,Beta,N,M1,M2)

    e_dof = nodal_dof * 2;
    Ke = zeros(e_dof,e_dof,num_e);
    A_I = A/I;

    for n = 1:num_e

        EI_Lo = (E*I)/Le_o(n);

        mat1 = [A_I, 0, 0; ...
                 0,  4, 2; ...
                 0,  2, 4];

        D = EI_Lo .* mat1;

        Le_w = Le(n);
        c = cos(Beta(n));
        s = sin(Beta(n));

        B_w = B(:,:,n);
        N_w = N(n);
        M1_w = M1(n);
        M2_w = M2(n);

        z = [s, -c, 0, -s, c, 0]';
        r = [-c, -s, 0, c, s, 0]';

        Ke(:,:,n) = B_w' * D * B_w ...
                  + (N_w/Le_w) * (z*z') ...
                  + (M1_w + M2_w)/(Le_w^2) * (r*z' + z*r');

    end

end

function [Fint] = global_internal_forcing(num_e,nodal_dof,B,N,M1,M2)

    num_nodes = num_e + 1;
    num_dof = num_nodes * nodal_dof;
    Fint = zeros(num_dof,1);

    for n = 1:num_e

        B_w = B(:,:,n);

        fVec = [N(n); M1(n); M2(n)];
        Fe_int = B_w' * fVec;

        node1 = n;
        node2 = n+1;
        dof_map = [3*node1-2, 3*node1-1, 3*node1, ...
                   3*node2-2, 3*node2-1, 3*node2];

        Fint(dof_map) = Fint(dof_map) + Fe_int;

    end

end

function [Kk] = formulate_global_stiffness_matrix(num_e,num_dof,Ke)

    Kk = zeros(num_dof,num_dof);

    for n = 1:num_e

        % Element rotation is accounted for by the co-rotational formulation.
        Ke_w = Ke(:,:,n);

        node1 = n;
        node2 = n+1;
        dof_map = [3*node1-2, 3*node1-1, 3*node1, ...
                   3*node2-2, 3*node2-1, 3*node2];

        Kk(dof_map,dof_map) = Kk(dof_map,dof_map) + Ke_w;

    end

end

function [X,Y,theta] = discretise_arch(Beta_A,c_A,L_beam,num_points)

    % Models an arch using the formulation from L.N. Virgin.
    % Geometry is then discretised with the specified number of points.

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
