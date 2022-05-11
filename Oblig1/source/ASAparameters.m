function P = ASAparameters;
%
% Defines parameters for the ASA
%

% Simulation grid
Nx = 1024;       % No of pixels, x-direction
Ny = 1024;       % No of pixels, y-direction

% Medium / problem
fc = 3e6;       % Center frequency, probe
c = 1500;       % Speed of sound
lambda = c/fc;  % Wavelength

% Physical grid
dx = lambda/2;  % pixel size, x-direction, [m]
dy = lambda/2;  % pixel size, y-direction, [m]
Dx = Nx*dx;     % Physical length, x-direction, [m]
Dy = Ny*dy;     % Physical length, y-direction, [m]

% Frequency grid
u = (2*(0:Nx-1)/Nx-1); % Defines as in fft
v = (2*(0:Ny-1)/Ny-1); % Aliasing will occure if sampling distance too large

% Axis
ax = linspace(-Dx/2, Dx/2, Nx);
ay = linspace(-Dy/2, Dy/2, Ny);

%% %%%%%%%%%
% Putting everything into a struct to be returned
P.Nx = Nx;
P.Ny = Ny;
P.fc = fc;
P.c  = c;
P.lambda = lambda;
P.dx = dx;
P.dy = dy;
P.Dx = Dx;
P.Dy = Dy;
P.u = u;
P.v = v;
P.ax = ax;
P.ay = ay;
