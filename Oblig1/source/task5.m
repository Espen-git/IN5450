clear all
diary on

P = ASAparameters;
lmd = P.lambda;
M = 24;
L = 1001;
theta = linspace(-pi/2,pi/2,L);
u = sin(theta);
kx = 2*pi/lmd*u;
d = lmd/2;
xpos = (1:M)*d;

%%
betas = 0:0.5:5;
k = kaiser(M, betas(1));
k = k/sum(k); % scale

W = beampattern(xpos, kx, k);
ana = analyzeBP(u, W);
three_db = ana.Three_dB; six_db = ana.Six_dB; max_sl = ana.maxSL; white_noise = 1/norm(k);

for i = 2:11
    k = kaiser(M, betas(i));
    k = k/sum(k); % scale
    W = beampattern(xpos, kx, k);
    ana = analyzeBP(u, W);
    three_db(i) = ana.Three_dB; six_db(i) = ana.Six_dB; max_sl(i) = ana.maxSL; white_noise(i) = 1/norm(k);
end

%%
subplot(2,2,1)
plot(betas, three_db);
xlabel('Beta')
ylabel('Width (theta)')
title('-3 dB Main lobe width')

subplot(2,2,2)
plot(betas, six_db);
xlabel('Beta')
ylabel('Width (theta)')
title('-6 dB Main lobe width')

subplot(2,2,3)
plot(betas, max_sl);
xlabel('Beta')
ylabel('Hight (dB)')
title('Maximum side lobe level')

subplot(2,2,4)
plot(betas, db(abs(white_noise)));
xlabel('Beta')
ylabel('WN Gain')
title('Kaiser window white noise gain')
