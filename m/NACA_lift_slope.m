function [a_0,alpha_L0] = NACA_lift_slope(NACA_num,c,N)
% This function finds the lift slope and the zero-lift angle of attack
% for a given airfoil with a specified chord length and N number of
% panels
%
% Inputs:   
%           NACA_num    - NACA airfoil number
%           c           - chord length
%           N           - number of employed panels to model the airfoil
%
% Outputs:  
%           a_0         - lift slope
%           alpha_L0    - zero-lift angle of attack
%
% Author: Cole MacPherson
% Collaborators: R. Block, Z. Lesan, S. Mansfield, A. Uprety
% Date: 27th Mar 2021

    %% Define m, p, and t of airfoil
    m = str2double(NACA_num(1))/100; % maximum camber
    p = str2double(NACA_num(2))/10; % location of maximum camber
    t = str2double(NACA_num(3:4))/100; % thickness

    %% Define airfoil points
    [x,y] = NACA_Airfoils(m,p,t,c,N); % x and y coordinates of 
    
    %% Solve for the coefficient of lift at varying AoA
    AoA = deg2rad(-5:10); % Define angles of attack
    c_l_vec = zeros(length(AoA),1); % initialize sectional coefficient of lift vector
    for i = 1:length(AoA)
        [c_l_vec(i),~] = Vortex_Panel(x,y,AoA(i)); % find sectional coefficient of lift
    end
    
    %% Define linear model of lfit slope
    a = polyfit(AoA,c_l_vec,1); % linear equation coefficients
    a_0 = a(1); % lift slope [1/rad]
    alpha_L0 = -a(2)/a(1); % zero-lift angle of attack [rad]
    
end