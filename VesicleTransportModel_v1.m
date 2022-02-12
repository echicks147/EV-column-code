function [output] = VesicleTransportModel_v1(param)
% =========================================================================
% Vesicle Transport Model
% =========================================================================
% Nicholas M.K. Rogers and Ethan C. Hicks
% Duke University, Durham, N.C.
% Last Updated: February 2022
% =========================================================================
% Scope: To model the distribution of extracellular vescicles (EVs) through sediment.
% =========================================================================
% Instrucitons: 
% =========================================================================
clear % to clear-out variable still in Workspace between runs

% -------------------------------------------------------------------------
% Input parameters (User action required):
% -------------------------------------------------------------------------
param = [3.65*10^-4, 100*10^-9, 25, 1000, 1135, 10^-3, 1.7*10^-4, 0.32, 5*10^-20, 1.03*10^21, 0, .00023];
dc        = param(1);   % Collector diameter (m)
dp        = param(2);   % Particle diameter (m)
T         = param(3);   % temp. (oC)
rho       = param(4);   % medium density (kg/m/m/m)
rho_p     = param(5);   % particle density (kg/m/m/m)
mu        = param(6);   % Dynamic or absolute viscosity (Pa*s)
v0        = param(7);   % Superficial velocity (m/s)
epsilon   = param(8);   % Porosity
A         = param(9);   % Hamaker constant (J)
n(1)      = param(10);  % Initial number concentration of EVs (#/m/m/m)
z(1)      = param(11);  % Initial distance (taken to be 0 as "at surface")
alpha     = param(12);  % Attachment efficiency (#/#)

% -------------------------------------------------------------------------
% Euler Integration (User action required):
% -------------------------------------------------------------------------
% Distance Parameter
startingz = -1; % Starting position (m) (usually 0 as a reference point)
endingz = 0; % Ending position (m) (depth, with surface as reference point)
distance = endingz - startingz; % Total distance calculation (m)
% Euler Integration steps
M = 10000; % Number of steps
zstep = distance/M; % height/depth of each step

% -------------------------------------------------------------------------
% Model Selection (User action required):
% -------------------------------------------------------------------------
% Select model for calculation by subsituting any positive intiger in place
% of 0
Levich = 0;
Yao = 0;
Rajagopalan_Tien = 0;
Tufenkji_Elimelech = 1;

% -------------------------------------------------------------------------
% Constants:
% -------------------------------------------------------------------------
g = 9.81; % gravitation constant (m/s/s)
kB = 1.38 * 10^(-23); % Boltzmann constant (m*m*kg/s/s/K)
T_K = T+273.15;

% -------------------------------------------------------------------------
% Dimensionless Numbers
% -------------------------------------------------------------------------
N_R = dp/dc; % Aspect Ratio
N_Lo = 4*A/(9*pi*dp^2*v0); % London Attraction Number
N_G = dp^2*(rho_p-rho)*g/(18*mu*v0); % Gravitational Number
N_Pe = 3*mu*dp*v0*dc/(kB*T_K); % Peclet Number
gamma = (1-epsilon)^(1/3); % Filter Coefficient
A_s = 2*(1-gamma^5)/(2-3*gamma+3*(gamma^5)-2*(gamma^6)); %Porosity-dependent Number
N_vdw = A/(kB*T_K); % van der Waals Number
% N_A = A/(12*pi*mu*(dp/2)*(dp/2)*v0); % Attraction Number

% -------------------------------------------------------------------------
% Eta calculations 
% -------------------------------------------------------------------------
% Levich (1962)
if Levich  > 0
Eta_B = 0.9*(kB*T_K/(mu*dp*dc*v0))^(2/3); % Brownian motion term
Eta_I = 0; % Shear motion term
Eta_G = 0; % Gravitational term
end
% Yao (1969)
if Yao > 0
Eta_B = 0; % Brownian motion term
Eta_I = (3/2)*((dp/2)/(dc/2))^2; % Shear motion term
Eta_G  = (rho_p-rho)*G*dp*dp/(18*mu*v0); % Gravitational term
end
% Rajagopalan and Tien (1976)
if Rajagopalan_Tien > 0
Eta_B = 4.04*A_s^(1/3)*N_Pe^(-2/3); % Brownian motion term
Eta_I = 0.72*A_s*(N_Lo^(1/8))*(N_R^(15/8)); % Shear motion term
Eta_G = (2.4*10^-3)*A_s*(N_G^1.2)*(N_R^-0.4); % Gravitaional term
end
% Tufenkji & Elimelech (2004)
if Tufenkji_Elimelech > 0
Eta_B = 2.4*(A_s^(1/3))*(N_R^-0.081)*(N_Pe^-0.715)*(N_vdw^0.052); % Brownian motion term
Eta_I = 0.55*A_s*(N_R^1.55)*(N_Pe^-0.125)*(N_vdw^0.125); % Shear motion term
Eta_G = 0.22*(N_R^-0.241)*(N_G^1.11)*(N_vdw^0.053); % Gravitational term
end
% Eta R calculation
Eta_T = Eta_B + Eta_I + Eta_G;
Eta_R = alpha*Eta_T;

% -------------------------------------------------------------------------
% Population Balance(s)
% -------------------------------------------------------------------------
for m = 1:M
% Linearized variables a and b
a = (-3/2)*((1-epsilon)*Eta_R/dc);
% Calculating change in population of free EVs as a function of depth, or
% distance (z).
z(m+1) = z(m)+zstep;
n(m+1) = n(m)+zstep*(a.*n(m));
end
[output] = (n);

% -------------------------------------------------------------------------
% Plots
% -------------------------------------------------------------------------
hold on
plot(z,n/n(1),'r','LineWidth', 1.5)
set(gcf,'color','w')
set(gca,'fontsize',20)
set(gca,'TickDir','out')
ylabel("n/n0")
xlabel("Distance (m)") 
title('1/20 PF-EV pH7 Crude Prep. no HA \alpha=0.00023')
xlim([0 1])
ylim([0 1])
%end
end
% =========================================================================
% END OF CODE
% =========================================================================