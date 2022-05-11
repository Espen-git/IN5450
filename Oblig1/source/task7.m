clear all
diary on

P = ASAparameters;
lmd = P.lambda;
M = 24;
L = 1001;
weights = ones(1,M)*(1/M);
theta = linspace(-pi/2,pi/2,L);
u = sin(theta);
kx = 2*pi/lmd*(u-(sin(3/4*pi)))*lmd;

n = 1:2:24;
d = 1/2;
e_n = [-0.017 -0.538 -0.617 -1.0 -1.142 -1.372 -1.487 -1.555 -1.537 -1.3 -0.772 -0.242];
d_n = (n/2 + e_n)*d;
ElPos = [-fliplr(d_n) d_n];

W = beampattern(ElPos, kx, weights);

subplot(2,1,1);
plot(u, db(abs(W)));
xlabel('u (sin(theta))')
ylabel('Response (db)')
title('3/4*pi Stearing')

% Unsteared -2 to 2
u2 = sin(theta)*2;
kx2 = 2*pi/lmd*(u)*lmd;
W2 = beampattern(ElPos, kx2, weights);

subplot(2,1,2);
plot(u2, db(abs(W2)));
xlabel('theta')
ylabel('Response (db)')
title('Unsteared for u in [-2,2]')
