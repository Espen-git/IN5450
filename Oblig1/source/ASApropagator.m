function Prop = ASApropagator(P,z);

% Propogator in ASA, from z_0 to z

Prop = zeros(P.Nx,P.Ny);

uSq = repmat(((P.u).^2)',1,P.Ny);
vSq = repmat((P.v).^2,P.Nx,1);

Prop = exp(1j*2*pi*z*sqrt(1-uSq - vSq)/P.lambda);

Prop(uSq+vSq>1) = 0;

