clear

%% Fixed stuff
L0 = 1.901e-3;       % Length scale [m]
t0 = 7.28e6^-1;      % Time scale [s]
e = 1.60217662e-19;  % Electron charge
kB = 1.38064852e-23; % Boltzmann constant

%% Simulation parameters direct computation
T0 =50;                 % Temperature (eV)
mi = 3.35e-27;       % Ionic mass (kg)
me = 9.109e-31;    % electronic mass (kg)
B0 = 1;                   % Magnetic field at magnetic axis (Tesla)
D = 5;                 % Diffusion in m^2/s
k0 = 2000;             % Stangeby
n0 = 1e19;             % Reference density
eps0 = 8.85e-12;


%% Simulation parameters inverse computation
a_inv = 2;           % Coefficient of the momentum equation
Drt_inv = 1;         % Coefficient of the drift velocity

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

% Reference speed
u0 = L0/t0;

% Compressibility coefficient
Mref = T0*e/(mi*u0^2);

% k*
k_star = t0^3*T0^(7/2)/(n0*L0^4)*k0/mi;

% tau_ie
% tau_ie = mi/me*2.4/3*1e10*T0^(0.5)*u0^2/n0/t0;
coef = 3*sqrt(2)/e^4/12*eps0^2*pi^1.5*mi/me*sqrt(me)*e^1.5;
tau_ie = 2/3*coef*T0^(0.5)*u0^2/n0/t0*mi/kB;

tau_ie*kB/mi*T0^1.5/n0

%% Results inverse computation

% Sound speed (m/s) 
cs_inv = sqrt(a_inv*(L0/t0)^2);

% Background temperature (eV) 
T0_inv = cs_inv^2/2/e*mi;

% Magnetic field at magnetic axis (Tesla)
B0_inv = 2*T0_inv*t0/L0^2/Drt_inv;
