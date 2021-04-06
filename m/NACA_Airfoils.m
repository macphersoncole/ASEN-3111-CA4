function [x,y] = NACA_Airfoils(m,p,t,c,N)
% Function to find the x and y coordinates for the surface of a NACA
% airfoil, whose contour is approzimated by N panels
%
% Inputs:   
%           m - maximum camber for the airfoil
%           p - location of the maximum camber of the airfoil
%           t - thickness of the airfoil
%           c - chord length of the airfoil
%           N - number of employed panels to model the airfoil
%
% Outputs:  
%           x - vector containing the x-location of boundary points
%           y - vector containing the y-location of boundary points
%
% Author: Cole MacPherson
% Collaborators: R. Block, Z. Lesan, S. Mansfield, A. Uprety
% Date: 27th Mar 2021

    %% Create x vector for x-values of each panel
    x_axis = 0:(c/(N/2)):c;

    %% Calculate the thickness distribution of the airfoil from the mean camber line
    y_t = ((t/0.2)*c)*...
        (0.2969.*sqrt(x_axis./c)...
        -0.126.*(x_axis./c)...
        -0.3516.*(x_axis./c).^2 ...
        +0.2843.*(x_axis./c).^3 ...
        -0.1036.*(x_axis./c).^4);
    
    %% Calculate the mean camber line
    if m ~= 0 && p ~= 0 % check if m and p are zero
        y_c = zeros(1,length(x_axis));
        for i = 1:length(x_axis) % interate over the whole of x_axis
            if x_axis(i) <= p*c % if the x position is before the maximum camber point
                y_c(i) = m*(x_axis(i)/(p^2))*(2*p-(x_axis(i)/c));
            elseif x_axis(i) > p*c % if the x position is after the maximum camber point
                y_c(i) = m*((c-x_axis(i))/((1-p)^2))*(1+(x_axis(i)/c)-2*p);
            end
        end
    else
        y_c = zeros(1,length(x_axis));
    end
    
    %% Calculate zeta
    zeta = atan(diff(y_c));
    zeta = [zeta, 0];
    
    %% Calculate the upper surface coordinates
    x_u = x_axis - y_t.*sin(zeta);
    y_u = y_c + y_t.*cos(zeta);
    
    %% Calculate the lower surface coordinates
    x_l = x_axis + y_t.*sin(zeta);
    y_l = y_c - y_t.*cos(zeta);
    
    %% Combine the upper and lower surfaces
    x = [flip(x_l), x_u(2:end)];
    y = [flip(y_l), y_u(2:end)];
    
end