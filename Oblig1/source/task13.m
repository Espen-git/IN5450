% Same as 12.m but with steering
clear all
diary on

P = ASAparameters;
lmd = 1;
M = 24;
L = 1001;
weights = ones(1,M)*(1/M);
theta = linspace(-pi/2,pi/2,L);
steering_angle = pi/3;
u = sin(theta);
kx = 2*pi/lmd*u;
kx_s = 2*pi/lmd*(u - (sin(steering_angle)));

%%
d = lmd;
xpos = (1:M)*d;
W1_e = (sin(kx.*d./2))./(kx./2);
W1 = beampattern(xpos, kx_s, weights);
W1_tot = W1_e' .* W1;
Result3 = analyzeBP(u,W1_tot);

%%
d = 2*lmd;
xpos = (1:M)*d;
W2_e = (sin(kx.*d./2))./(kx./2);
W2 = beampattern(xpos, kx_s, weights);
W2_tot = W2_e' .* W2;
Result4 = analyzeBP(u,W2_tot);

%%
figure();
subplot(1,2,1);
plot(u, db(abs(W1_tot)));
ylim([-50 10]);
xlabel('u [sin(theta)]');
ylabel('Response [dB]');
title('Beampattern for d=lmd')

subplot(1,2,2);
plot(u, db(abs(W2_tot)));
ylim([-50 10]);
xlabel('u [sin(theta)]');
ylabel('Response [dB]');
title('Beampattern for d=2*lmd')

%%
% Element factors
figure();
subplot(1,2,1);
plot(u, db(abs(W1_e)));
ylim([-50 10]);
xlabel('u [sin(theta)]');
ylabel('Response [dB]');
title('EF for d=lmd')

subplot(1,2,2);
plot(u, db(abs(W2_e)));
ylim([-50 10]);
xlabel('u [sin(theta)]');
ylabel('Response [dB]');
title('EF for d=2*lmd')