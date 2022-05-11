function U0 = ASAsource(P,Type,varargin);

% Function setting up the source for the ASA simulation

if strcmp(Type,'piston')
    % Defining a piston w/uniform pressure, radius 'r' placed in origo
    r = varargin{1};
    [x,y] = meshgrid(P.ax, P.ay);
    U0 = zeros(P.Nx,P.Ny);
    U0(x.^2 + y.^2 <= r^2) = 1;
elseif strcmp(Type,'square')
    % Defining a piston w/uniform pressure, radius 'r' placed in origo
    r = varargin{1};
    [x,y] = meshgrid(P.ax, P.ay);
    U0 = zeros(P.Nx,P.Ny);
    U0(abs(x) <= r & abs(y) <= r) = 1;
else
    error('Source type not defined in ASAsource.m')
end



