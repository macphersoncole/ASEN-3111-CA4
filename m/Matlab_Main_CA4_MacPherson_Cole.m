%% ASEN 3111 - Computational Assignment 4 - Main
% Script to model the flow over a finite wing using Prandtl Lifting Line
% Theory (PLLT) to calculate lift and drag. Additionally, an error analysis
% study was conducted to determine the affect of the number of panels used 
% in the PLLT function. Lastly, see how design facotrs like taper and
% aspect ratio affect incuded drag and aerodynamic efficiency
%
% Author: Cole MacPherson
% Collaborators: R. Block, Z. Lesan, S. Mansfield, A. Uprety
% Date: 27th Mar 2021

%% Housekeeping

clc;
clear;
close all;
tic

%% Problem 2 (e, C_L, and, C_Di for specified wing & error analysis of PLLT)
% define variables
b = 100; % span [ft]
c_r = 15; % root chord [ft]
c_t = 5; % tip chord [ft]
NACA_r = '2412'; % root airfoil
NACA_t = '0012'; % tip airfoil
geo_r = deg2rad(5); % root geometric angle of attack [rad]
geo_t = 0; % root geometric angle of attack [rad]
V_inf = 220; % free-stream velocity [ft/s]
rho_inf = 0.0023769; % free-stream air density [slugs/ft^3]
S = (1/2)*(c_t+c_r)*b; % surface area [ft^2]

% Calculate lift slope and zero-lift angle of attack of each NACA airfoil
[a0_r, aero_r] = NACA_lift_slope(NACA_r,c_r,100); % NACA root airfoil
[a0_t, aero_t] = NACA_lift_slope(NACA_t,c_t,100); % NACA tip airfoil

% Calculate e, c_L, and C_di
[e,C_L,C_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,5000);
fprintf('Problem 2:\n\te:    %0.4f\n\tc_L:  %0.4f\n\tc_Di: %0.4f\n\n',e,C_L,C_Di);

% N required for less than 5%, 1%, and 0.1% error
    % initialization
    N = 2; 
    i = 1;
    err = 100;
    % while loop to find the error and number of odd terms associated with that
    % error
    N_vec = zeros(32,1); % Preallocated for speed because I know how large in need to be based on previous runs
    c_L_vec = zeros(32,1);
    c_D_vec = zeros(32,1);
    err_vec = zeros(32,1);
    while err > 0.05
        N_vec(i) = N;
        [~, c_L_vec(i), c_D_vec(i)] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N);
        err_vec(i) = 100*abs(C_L - c_L_vec(i))/C_L;
        err = err_vec(i);
        i = i + 1;
        N = N + 1;
    end
    % for loop to find the N value associated with desired error values
    for err = [5, 1, 0.1]
      index = find(err_vec <= err, 1);
      L = (1/2)*c_L_vec(index)*rho_inf*S*V_inf^2;
      D = (1/2)*c_D_vec(index)*rho_inf*S*V_inf^2;
      fprintf('\tFor a less thsn %.1f%% error: %i odd terms are required\n\t\t(L = %0.4f lbf)\n\t\t(D = %0.4f lbf)\n',err,N_vec(index),L,D);
    end
    fprintf('\n');

    % plot error vs number of odd terms
    figure
    plot(N_vec,err_vec)
    grid on
    title('Number of Odd Terms vs Percent Error');
    xlabel('Number of Odd Terms');
    ylabel('Percent Error');

%% Problem 3 (design factors affect on induced drag and aerodynamic efficiency)
AR = 4:2:10; % AR values
taper_ratios = 0:0.01:1; % taper ratio values
e_vals = zeros(length(AR),length(taper_ratios)); % initialize matrix that will store span efficiency values
k = 1;
% find the span efficiency values for each AR and taper ratio
for i = AR
    l = 1;
    for j = taper_ratios
        c_t = c_r*j;
        b = i*((1+j)*c_r)/2;
        [e_vals(k,l),~,~] = PLLT(b,2*pi,2*pi,c_t,c_r,0,0,1,1,N);
        l = l + 1;
    end
    k = k + 1;
end

% plot results
figure
hold on
for i = 1:4
    plot(taper_ratios,e_vals(i,:),'DisplayName',sprintf('AR = %i',AR(i)));
end
title('Span Efficiency Factor vs Taper Ratio');
xlabel('Taper Ratio [c_t/c_r]');
ylabel('e');
legend('location','southeast');

figure
hold on
for i = 1:4
    plot(taper_ratios,(1./e_vals(i,:))-1,'DisplayName',sprintf('AR = %i',AR(i)));
end
title('Induced Drag Factor vs Taper Ratio');
xlabel('Taper Ratio [c_t/c_r]');
ylabel('\delta');
legend('location','northeast');

%% End Housekeeping
toc
