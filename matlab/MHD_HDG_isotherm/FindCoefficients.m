clear

%% Fixed stuff
L0 = 1.901e-3;       % Length scale
t0 = 7.28e6^-1;      % Time scale
e = 1.60217662e-19;  % Electron charge

%% Simulation parameters direct computation
T0 =50;             % Temperature (eV)
mi = 3.35e-27;  % Ionic mass (kg)
B0 = 1;              % Magnetic field at magnetic axis (Tesla)
D = 1;                % Diffusion in m^2/s


%% Simulation parameters inverse computation
a_inv = 25;           % Coefficient of the momentum equation
Drt_inv = 1;         % Coefficient of the drift velocity
D_ad_inv = 0.003;        % Adimensional diffusion

%% Results direct computation
% Plasma sound speed (m/s)
cs = sqrt(2*T0*e/mi);

% Coefficient of the momentum equation
a = cs^2/(L0/t0)^2; 

% Adimensional diffusion
Dstar = D*t0/L0^2;

% Coefficient of the drift velocity
Drt = 2*T0*t0/(L0^2*B0);

% Larmor rayon 
rho_l = sqrt(mi*T0/(e*B0^2));

% Ciclotronic frequency
omega = e*B0/mi;

D_ad = D*t0/L0^2;

%% Results inverse computation

% Sound speed (m/s) 
cs_inv = sqrt(a_inv*(L0/t0)^2);

% Background temperature (eV) 
T0_inv = cs_inv^2/2/e*mi;

% Magnetic field at magnetic axis (Tesla)
B0_inv = 2*T0_inv*t0/L0^2/Drt_inv;

% Physical diffusion
D_phy_inv = D_ad_inv *L0^2/t0;
