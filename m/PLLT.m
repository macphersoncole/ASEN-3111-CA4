function [e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)
% This function solves the findemental equation of Prandtl Lifting Line
% Theory for fintie wings with thick airfoils, whose contour is 
% approzimated by N votex panels of linearly varying strength.
%
% Inputs:
%           b       - span
%           a0_t    - cross-sectional lift slope at the tip [1/rad]
%           a0_r    - cross-sectional lift slope at the root [1/rad]
%           c_t     - chord length at the tip
%           c_r     - chord length at the root
%           aero_t  - zero-lift angle of attack at the tip [rad]
%           aero_r  - zero-lift angle of attack at the root [rad]
%           geo_t   - geometric angle of attack at the tip [rad]
%           geo_r   - geometric angle of attack at the root [rad]
%           N       - number of odd terms to include in the seies 
%                       expansion for circulation
% Outputs:
%           e       - span efficiency factor
%           c_L     - coefficient of lift
%           c_Di    - coefficient of induced drag
%
% Author: Cole MacPherson
% Collaborators: R. Block, Z. Lesan, S. Mansfield, A. Uprety
% Date: 27th Mar 2021
    
    %% Calculate Wing Aspects
    S = (1/2)*(c_t+c_r)*b; % surface area of the wing [ft^2]
    AR = (b^2)/S; % aspect ratio of the wing
    theta = ((1:N)'.*pi)./(2*N); % theta values as defined in lab doc (half of the wing) [rad]
    c = c_r + (c_t-c_r).*cos(theta); % linear chord length variation from tip to root [ft]
    a0 = a0_r + (a0_t-a0_r).*cos(theta); % linear cross-sectional lift slope variation from tip to root [1/rad]
    alpha_L0 = aero_r + (aero_t-aero_r).*cos(theta); % linear zero-lift angle of attack variation from tip to root [rad]
    alpha_geo = geo_r + (geo_t-geo_r).*cos(theta); % linear geometric angle of attack variation from tip to root [rad]
    n = (1:2:2*N)'; % initialize n values
    
    %% Calculate Angle of Attack
    alpha = alpha_geo - alpha_L0; % [rad]
    
    %% Calculate the A_n Coefficients
    C = zeros(N,N);
    for i = 1:N
        for j = 1:N
            C(i,j) = 4*b*sin(n(j)*theta(i))/(a0(i)*c(i))...
                + n(j)*(sin(n(j)*theta(i))/sin(theta(i))); % A_n coefficient matrix (A_n coefficient at each theta value)
        end
    end
    
    %% Calculate A_n terms
    A_n = C\alpha; % A_n values
    
    %% Calculate delta
    delta = sum(n(2:end).*((A_n(2:end)./A_n(1)).^2)); % delta
    
    %% Calculate Span Efficiency Factor
    e = (1+delta)^(-1); % span efficiency factor
    
    %% Calculate Coeficient of Lift
    c_L = A_n(1)*pi*AR; % coefficient of lift
    
    %% Calculate Coeficient of Induced Drag
    c_Di = (c_L^2)/(pi*e*AR); % coefficient of induced drag

end