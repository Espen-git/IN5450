clear all
diary on

L = 101;
angels = linspace(-pi/3, pi/3, L);
anas = [];

for angle = angels
P = ASAparameters;
lmd = P.lambda;
M = 24;
weights = ones(1,M)*(1/M);
theta = linspace(-pi/2,pi/2,5001);
u = sin(theta);
kx = 2*pi/lmd*(u-sin(angle));

n = 1:2:24;
d = 1/2*lmd;
e_n = [-0.017 -0.538 -0.617 -1.0 -1.142 -1.372 -1.487 -1.555 -1.537 -1.3 -0.772 -0.242];
d_n = (n/2 + e_n)*d;
ElPos = [-fliplr(d_n) d_n];

W = beampattern(ElPos, kx, weights);
ana = analyzeBP(u, W);
anas = [anas, ana];
end
%
values = struct2cell(anas);
values = cell2mat(values);
tree_db = reshape(values(1,:,:),[1,L]);
six_db = reshape(values(2,:,:),[1,L]);
subplot(1,2,1);
hold on
plot(angels, sin(tree_db));
plot(angels, sin(six_db));
xlabel('angle'); ylabel('Width [sin(angle)]');
title('-3 dB, -6 dB as function of steering angle');
legend('-3 dB', '-6 dB');

subplot(1,2,2);
hold on
plot(angels, tree_db);
plot(angels, six_db);
xlabel('angle'); ylabel('Width [angle]');
title('-3 dB, -6 dB as function of steering angle');
legend('-3 dB', '-6 dB');
