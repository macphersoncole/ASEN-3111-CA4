function [c_l, C_p] = Vortex_Panel(x,y,alpha)
% Function to compute the sectional coefficient of lift and coefficient of
% pressure around an airfoil, whose contour is approzimated by N vortex
% panels of linearly varying strength
%
% Inputs:   
%           x       - x-coordinates of airfoil
%           y       - y-coordinates of airfoil
%           V_inf   - free-stream velocity of the flow
%           alpha   - angle of attack of the airfoil
%
% Outputs:  
%           c_l     - sectional coefficient of lift for the airfoil
%           C_p     - coefficient of pressure for the airfoil
%
% Author: Cole MacPherson
% Collaborators: R. Block, Z. Lesan, S. Mansfield, A. Uprety
% Date: 27th Feb 2021
    
    %% Declare constants
    rho = 1.225; % density of air at sea level [kg/m^3]
    c = max(x) - min(x); % chord length [m]
    
    %% Variable initialization
    N = length(x) - 1; % number of panels
    X = zeros(N,1); % will store the x-length of each panels
    Y = zeros(N,1); % will store the y-length of each panles
    S = zeros(N,1); % will store the length of each panel
    theta = zeros(N,1); % will store the theta value for each panel
    sine = zeros(N,1); % will store the sine value evaluated at theta
    cosine = zeros(N,1); % will store the cosine values evaluated at theta
    RHS = zeros(N,1); % will store the RHS of system of equation
    CN1 = zeros(N,N); % will store normal coefficient 1
    CN2 = zeros(N,N); % will store normal coefficient 2
    CT1 = zeros(N,N); % will store tangential coefficient 1
    CT2 = zeros(N,N); % will store tangential coefficient 2
    AN = zeros(N+1,N+1); % will store the normal-velocity influence coefficient
    AT = zeros(N,N+1); % will store the tangential-velocity influence coefficient
    ds = zeros(N,1); % will store the infinitesimal point
    V = zeros(N,1); % will store the velocity
    C_p = zeros(N,1); % will store the coefficient of pressure
    
    %% Compute coordinates (X,Y) of each control point and panel length S for each panel
    for i = 1:N
        X(i) = (1/2)*(x(i)+x(i+1)); % x-location of panel
        Y(i) = (1/2)*(y(i)+y(i+1)); % y-loaction of panel
        S(i) = sqrt((x(i+1)-x(i))^2 + (y(i+1)-y(i))^2); % length of panel
        theta(i) = atan2(y(i+1)-y(i),x(i+1)-x(i)); % theta value at each panel
        sine(i) = sin(theta(i)); % sine value corresponding to the panels theta value
        cosine(i) = cos(theta(i)); % cosine value corresponding to the panels theta value
        RHS(i) = sin(theta(i)-alpha); % right hand side of the equation
    end
    
    %% Compute tangential and normal coefficients
    for i = 1:N
        for j = 1:N
            if i == j
                % tangential and normal coefficients if i equals j
                CN1(i,j) = -1;
                CN2(i,j) = 1;
                CT1(i,j) = pi/2;
                CT2(i,j) = pi/2;
            else
                % calculate constants
                A = -(X(i)-x(j))*cosine(j) - (Y(i)-y(j))*sine(j);
                B = (X(i)-x(j))^2 + (Y(i)-y(j))^2;
                C = sin(theta(i)-theta(j));
                D = cos(theta(i)-theta(j));
                E = (X(i)-x(j))*sine(j) - (Y(i)-y(j))*cosine(j);
                F = log(1+S(j)*(S(j)+2*A)/B);
                G = atan2(E*S(j),B+A*S(j));
                P = (X(i)-x(j))*sin(theta(i)-2*theta(j)) + ...
                    (Y(i)-y(j))*cos(theta(i)-2*theta(j));
                Q = (X(i)-x(j))*cos(theta(i)-2*theta(j)) - ...
                    (Y(i)-y(j))*sin(theta(i)-2*theta(j));
                % tangential and normal coefficients if i does not equal j
                CN2(i,j) = D + (1/2)*Q*(F/S(j)) - (A*C+D*E)*G/S(j);
                CN1(i,j) = (1/2)*D*F + C*G - CN2(i,j);
                CT2(i,j) = C + (1/2)*P*F/S(j) + (A*D-C*E)*G/S(j);
                CT1(i,j) = (1/2)*C*F - D*G - CT2(i,j);
            end
        end
    end
    
    %% Compute tangential and normal velocity influence coefficients
    for i = 1:N
        AN(i,1) = CN1(i,1);
        AN(i,N+1) = CN2(i,N);
        AT(i,1) = CT1(i,1);
        AT(i,N+1) = CT2(i,N);
        for j = 2:N
            AN(i,j) = CN1(i,j) + CN2(i,j-1);
            AT(i,j) = CT1(i,j) + CT2(i,j-1);
        end       
    end
    
    AN(N+1,1) = 1;
    AN(N+1,N+1) = 1;
    for j = 2:N
        AN(N+1,j) = 0;
    end
    RHS(N+1) = 0;
    
    %% Compute gamma
    gamma = AN\RHS; % vortex streangth
    
    %% Compute velocity and coefficient of pressure
    for i = 1:N
        V(i) = cos(theta(i)-alpha);
        for j = 1:N+1
           V(i) = V(i) + AT(i,j)*gamma(j); % velocity
           C_p(i) = 1 - (V(i))^2; % coefficient of pressure
        end
    end
    
    %% Calculate infinitesimal point
    for i = 1:N
        ds(i) = sqrt((x(i)-x(i+1))^2 + (y(i)-y(i+1))^2); % infinitesimal point
    end
    ds = [ds; 0];
    
    %% Calculate Gamma
    Gamma = sum(2*pi*gamma.*ds); % circulation
    
    %% Calculate sectional lift
    lift = rho*Gamma; % sectional lift
    
    %% Calculate sectional cefficient of lift
    c_l = lift/((1/2)*rho*c); % sectional coefficient of lift
    c_l = sum(c_l); % sum of the sectional coefficient of lift vector
    
end