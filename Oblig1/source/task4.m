clear all
diary on

P = ASAparameters;
lmd = P.lambda;
M = 24;
L = 1001;
weights = ones(1,M)*(1/M);
theta = linspace(-pi/2,pi/2,L);
u = sin(theta);
kx = 2*pi/lmd*u;

%%
d = lmd/4;
xpos = (1:M)*d;
W1 = beampattern(xpos, kx, weights);
Result1 = analyzeBP(u,W1)

%%
d = lmd/2;
xpos = (1:M)*d;
W2 = beampattern(xpos, kx, weights);
Result2 = analyzeBP(u,W2)

%%
d = lmd;
xpos = (1:M)*d;
W3 = beampattern(xpos, kx, weights);
Result3 = analyzeBP(u,W3)

%%
d = 2*lmd;
xpos = (1:M)*d;
W4 = beampattern(xpos, kx, weights);
Result4 = analyzeBP(u,W4)

%%
subplot(2,4,1);
plot(u, db(abs(W1)));
ylim([-30 0]);
xlabel('u [sin(theta)]');
ylabel('Response [db]');
title('Beampattern for d=lmd/4')

subplot(2,4,5);
plot(theta, db(abs(W1)));
ylim([-30 0]);
xlabel('[theta]');
ylabel('Response (db)');
title('Beampattern for d=lmd/4')

subplot(2,4,2);
plot(u, db(abs(W2)));
ylim([-30 0]);
xlabel('u [sin(theta)]');
ylabel('Response [db]');
title('Beampattern for d=lmd/2')

subplot(2,4,6);
plot(theta, db(abs(W2)));
ylim([-30 0]);
xlabel('[theta]');
ylabel('Response (db)');
title('Beampattern for d=lmd/2')

subplot(2,4,3);
plot(u, db(abs(W3)));
ylim([-30 0]);
xlabel('u [sin(theta)]');
ylabel('Response [db]');
title('Beampattern for d=lmd')

subplot(2,4,7);
plot(theta, db(abs(W3)));
ylim([-30 0]);
xlabel('[theta]');
ylabel('Response (db)');
title('Beampattern for d=lmd')

subplot(2,4,4);
plot(u, db(abs(W4)));
ylim([-30 0]);
xlabel('u [sin(theta)]');
ylabel('Response [db]');
title('Beampattern for d=2*lmd')

subplot(2,4,8);
plot(theta, db(abs(W4)));
ylim([-30 0]);
xlabel('[theta]');
ylabel('Response (db)');
title('Beampattern for d=2*lmd')
