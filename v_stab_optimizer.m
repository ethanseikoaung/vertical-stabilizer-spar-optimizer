clc, clear

% Spar Geometry
l = 0.49455; % Semi-spar length in [m]
h_range = 0.00635; % Spar height range in m
w_range = 0.001:0.0001:0.05; % Spar width range in m
cross_dimen = [];

% Material Properties
sig_max = 14e6; % Pa
tau_max = 1.5e6; % Pa
rho = 160; %kg/m^3

% Safety Factor
SF = 1.6; 

% Design Values
sig_des = sig_max/SF;
tau_des = tau_max/SF;


% Estimated Lift

Lmax = 20; % Maximum aerodynamic load in [N]

% Variable Initialization
L=[];
V=[];
M=[];
X = 0:0.001:l;

% Reactions at root of top spar
R_y = - Lmax*l*pi/4;
M_R_z = -Lmax*(l^2)/3;

fprintf("R_y: %6.5f [N]\n",R_y);
fprintf("M_R_z: %6.5f [Nm]\n",M_R_z);

for x = 0:0.001:l

    % Load distribution on the aft spar [N/m]
    L_x = Lmax*sqrt((1-(x./l).^2));
    L = [L,L_x];
    % Shear along the Spar
    V_x =  - integral(@(s) Lmax*sqrt((1-(s./l).^2)), x, l);
    V=[V,V_x];
    % Bending Moment about the z axis along the Spar
    f = @(s) (s-x).*Lmax.*sqrt((1-(s./l).^2));
    M_x = integral(f, x, l);
    M = [M,M_x];
end

%Spar Size Optimization
for i = 1:length(h_range)
    for j = 1:length(w_range)
        h = h_range(i);
        w = w_range(j);
        I = w*h^3/12; % Moment of Inertia about neutral axis
        Q = (w*h/2)*(h/4); % Moment of Area

        % Shear Stress at the neutral axis at the root
        tau =  abs(R_y)*Q/(I*w);

        % Bending stress at the outer surface
        sig = abs(M_R_z*h/(2*I));

        if tau < tau_des && sig < sig_des
            cross_dimen = [cross_dimen; [h,w,h*w]];
        end
    end
end

% Finding the dimension with the minimal cross-sectional area

if isempty(cross_dimen)
    error('No feasible spar geometry found. Relax the constraints or expand the search range.');
else
    [min_cross, index] = min(cross_dimen(:,3));
    opti_h = cross_dimen(index,1);
    opti_w = cross_dimen(index,2);
    min_weight = rho*l*opti_w*opti_h*9.81;
end

% Lift Distribution
figure
plot(X,L,'b-','LineWidth',1.5);
grid on;
xlabel('Spanwise Position, x (m)');
ylabel('Lift Distribution (N/m)');
title("Lift Distribution along Semi-span");

% Shear Diagram
figure 
plot(X,V,'b-','LineWidth',1.5);
grid on;
xlabel('Spanwise Position, x (m)');
ylabel('Shear Force (N)');
title("Shear Distribution along Aft Span");

% Bending Moment Diagram
figure 
plot(X,M,'r-','LineWidth',1.5);
grid on;
xlabel('Spanwise Position, x (m)');
ylabel('Bending Moment (N·m)');
title("Bending Moment along Aft Span");

%Optimal Geometry
fprintf('---------------Optimal Geometry---------------\n');
fprintf("Optimal Height: %6.6f [m] \n",opti_h);
fprintf("Optimal Width: %6.6f [m] \n",opti_w);
fprintf("Optimal Cross Sectional Area %6.6f [m^2] \n",min_cross);
fprintf("Weight: %6.6f [N] \n",min_weight);