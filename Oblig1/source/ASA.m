clear all
diary on
%
% Script that sets up and runs an anglular spectrum approach simulation
%

%%%%%%%%%%%%%%%%%%
% Parameters read from file

P = ASAparameters;

%%%%%%%%%%%%%%%%%%
% Source, i.e. U(x,y,z_0), defined in separate function

% U0 = ASAsource(P,'piston',15*P.lambda);

U0 = ASAsource(P,'square',15*P.lambda);

%%%%%%%%%%%%%%%%%%
% Define propagator for distance (z-z0)
z = 0.08; % [m]
Prop = ASApropagator(P,z);

U = ifft2t(fft2t(U0).*Prop);

%%
figure(); imagesc(P.ax,P.ay,db(abs(U))); 
CA = caxis; caxis(CA(2) + [-40 0]); colorbar;
xlabel('x [m]'); ylabel('y [m]'); 
title(['Wave field at depth z = ',num2str(100*z,'%.1f'),' cm. Estimated using ASA.']);

figure();
plot(P.ax,db(abs(U(:,floor(P.Nx/2)))),'LineWidth',2);
ylim(max(db(abs(U(:))))+[-40 0]); grid on;
xlabel('x [m]'); ylabel('[dB]');
title(['Center cut of wave field at depth z = ',num2str(100*z,'%.1f'),' cm. Estimated using ASA.']);

%% Smoothing with haning window
Kernel = hanning(3)*hanning(3).';
Kernel = Kernel/sum(Kernel(:));
U0Mod = conv2(U0,Kernel,'same');
UMod = ifft2t(fft2t(U0Mod).*Prop);

figure(); imagesc(P.ax,P.ay,db(abs(UMod))); 
CA = caxis; caxis(CA(2) + [-40 0]); colorbar;
xlabel('x [m]'); ylabel('y [m]'); 
title(['Wave field at depth z = ',num2str(100*z,'%.1f'),' cm. Estimated using ASA.']);

figure();
plot(P.ax,db(abs(UMod(:,floor(P.Nx/2)))),'LineWidth',2);
ylim(max(db(abs(UMod(:))))+[-40 0]); grid on;
xlabel('x [m]'); ylabel('[dB]');
title(['Center cut of wave field at depth z = ',num2str(100*z,'%.1f'),' cm. Estimated using ASA.']);



%% Task 2
z = 0.40;
dz = 2*P.lambda;
nZ = round(z/dz);
Zs = linspace(0,z,nZ);
Resp = zeros(P.Nx,nZ);
Prop = ASApropagator(P,dz);
U = ifft2t(fft2t(U0Mod).*Prop);
Resp(:,1) = U(:,floor(P.Ny/2)).';

for ii=2:nZ
    if mod(ii,25)==0
        fprintf('#');
    end
    
    U = ifft2t(fft2t(U).*Prop);
    Resp(:,ii) = U(:,floor(P.Ny/2)).';
end

%%
% 2D image of predicted wave field for y=0
figure(); 
imagesc(Zs, P.ax, db(abs(Resp)));
xlabel('z [m]'); ylabel('x [m]');
title('wave field for y=0')

% pressure of the center of field. Along z-axis.
figure(); 
plot(Zs, db(abs(Resp(floor(P.Nx/2), :))));
xlabel('z [m]'); ylabel('[dB]');
title('Pressure of the center of field. Along z-axis.')

%%
% Center cut near-field far-field limit.
d = ((15*P.lambda)^2)/(P.lambda);
[v, ix] = min(abs(Zs-d)); % ix -> index closest to near-field far-field limit.

figure();
plot(P.ax, db(abs(Resp(:,ix))),'LineWidth',2);
ylim(max(db(abs(UMod(:))))+[-40 0]); grid on;
xlabel('x [m]'); ylabel('[dB]');
title(['Center cut of wave field at depth z = ',num2str(d,'%.4f'),' cm. Estimated using ASA.']);

%% Center cut much further than the near-field far-field limit.
figure();
plot(P.ax, db(abs(Resp(:,end))),'LineWidth',2);
ylim(max(db(abs(UMod(:))))+[-40 0]); grid on;
xlabel('x [m]'); ylabel('[dB]');
title(['Center cut of wave field at depth z = ',num2str(Zs(end),'%.4f'),' cm. Estimated using ASA.']);


