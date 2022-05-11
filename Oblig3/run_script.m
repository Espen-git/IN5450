clear all
close all
load('mimo_project.mat')
filetype = '-depsc'; % eps with colour
load_data = 1; % 0 = load, 1 = calculate
dr = 70; % dynamic range of images
yres = 1000; % resolution of imagesclose all

xres = yres*2;

%% Plot resevers and transmitters positions. Plot and calculate viritual array positions.
figure('Position',[300 300 1200 350]);
hold on
rx_zeros = zeros(1, length(rx_pos));
tx_zeros = zeros(1, length(tx_pos));
plot(rx_pos, rx_zeros, 'k.')
plot(tx_pos(1), 0, 'bd')
plot(tx_pos(2), 0, 'rd')

% viritual array
% Split per transmitter only for plotting purposes
vx_t1_pos = (tx_pos(1) + rx_pos)/2;
vx_t2_pos = (tx_pos(2) + rx_pos)/2;
vx_pos = sort(cat(2, vx_t1_pos,vx_t2_pos));
vx_zeros = zeros(1, length(vx_pos));

plot(vx_t1_pos, vx_zeros(1:32), 'b|')
plot(vx_t2_pos, vx_zeros(1:32), 'r|')
legend('Resevers','Transmitter 1','Transmitter 2','Virtual elements from T1','Virtual elements from T2')
xlabel('X position [m]')
ylabel('Y position [m]')
grid('on')
hold off

print('images/Arrays',filetype)
% Slightly larger, but mailny twise the sampling rate.

% Theoretical latteral resolution = lambda / L = lambda m / (2*0.26775) m =
lambda = c/fc; % = 0.0340 m
% 0.0340 / 0.5355 = 0.0635 (radians) = 3.638282 (degrees)

% Approx lateral resolution at 4 m: L = R*a = 4 (m) * 0.0635 = 0.254 m
% Axial resolution = c / (2*B) = 0.0170 m
% At 4 meters 0.254 / 0.0170 = 14.9, axial is 15 times better than latteral.

%% Syntetic replicas and pulse compression
alpha = B/T_p; % chirp rate
t = [0:(1/fs):T_p]; % time for pulse
tmax = (N_t/fs - 1/fs); % end time for recever
tt =  [0:(1/fs):tmax]; % time for recever
LFM_up = exp(j * 2 * pi * ( (fc - (B/2)) * t + ((alpha * (t.^2))/2) ));
LFM_down = exp(j * 2 * pi * ( (fc + (B/2)) * t - ((alpha * (t.^2))/2) ));

tmp = tdma_data(:,1,1); % single receve, single transmit
[m,lag] = xcorr(tmp, LFM_up);
% only with positive lag
m = m(lag >= 0);
tm = lag(lag >= 0)/fs;  


%figure()
%plot(t,real(LFM_up))
%title('Transmitted Pulse (LFM upchirp)')
%ylabel('Magnitude (L.U.)')
%xlabel('[s]')
%print('images/transmitted_pulse',filetype)

figure()
ax2 = subplot(2,1,1);
plot(tt,abs(tmp),'LineWidth',2)
title('Received tdma (before pulse compression)')
xlabel('Time (s)')
ylabel('Magnitude (L.U.)')
ylim([0 0.025])

ax3 = subplot(2,1,2);
plot(tm,abs(m))
title('Received tdma (after pulse compression)')
xlabel('Time (s)')
ylabel('Magnitude (L.U.)')
ylim([0 48])

linkaxes([ax2,ax3],'x')
print('images/compression_tdma',filetype)

% CDMA compression
tmp = cdma_data(:,20); % single receve, both transmit
[m1,lag1] = xcorr(tmp, LFM_down);
[m2,lag2] = xcorr(tmp, LFM_up);
% only with positive lag
m1 = m1(lag1 >= 0);
m2 = m2(lag2 >= 0);
tm1 = lag1(lag1 >= 0)/fs;
tm2 = lag2(lag2 >= 0)/fs; 

figure()

subplot(3,1,1);
plot(tt,abs(tmp),'LineWidth',2)
title('Received cdma (before pulse compression)')
xlabel('Time (s)')
ylabel('Magnitude (L.U.)')
ylim([0 0.025])

subplot(3,1,2);
plot(tm1,abs(m1))
title('Received cdma (after pulse compression LFM\_down)')
xlabel('Time (s)')
ylabel('Magnitude (L.U.)')
ylim([0 48])

subplot(3,1,3);
plot(tm2,abs(m2))
title('Received cdma (after pulse compression LFM\_up)')
xlabel('Time (s)')
ylabel('Magnitude (L.U.)')
ylim([0 48])
print('images/compression_cdma',filetype)

% Theoretical time resolutoin = 1 / B = 1.0000e-04 = 0.0001 seconds
% Measured resolution = 0.01792 - 0.017725 = 1.9500e-04 = 0.000195 seconds
% Measured ~ 2 * Theoretical

%% DAS TDMA
% xcorr over all data
all_m = zeros(N_t, N_tx, N_rx);
%all_tm = zeros(N_t, N_tx, N_rx);

for transmit = 1:N_tx
    for receiver = 1:N_rx
        tmp_data = tdma_data(:,receiver,transmit); % single receve, single transmit
        [m,lag] = xcorr(tmp_data, LFM_up);
        m = m(lag >= 0);
        tm = lag(lag >= 0)/fs;
        all_m(:,transmit,receiver) = m;
        %all_tm(:,transmit,receiver) = tm;
    end
end

% DAS
x = linspace(-5,5,xres);
y = linspace(0,5,yres);
image_tdma_2_transmit = complex(zeros(yres,xres));
image_tdma_1_transmit = complex(zeros(yres,xres));

%load_data = 1; % 0 = load, 1 = calculate

if load_data == 0
    load('image_tdma_1_transmit.mat');
    load('image_tdma_2_transmit.mat');
elseif load_data == 1
    for transmit = 1:N_tx
        for receiver = 1:N_rx
            m = all_m(:,transmit,receiver);
            for xx = 1:xres
                for yy = 1:yres
                    r_t = sqrt( (tx_pos(transmit)-x(xx))^2 + (y(yy))^2 ); 
                    r_r = sqrt( (rx_pos(receiver)-x(xx))^2 + (y(yy))^2 );
                    r = r_t + r_r;
                    t_delay = r / c;
                    t_sample = round(t_delay * fs);
                    if t_sample > 0 && t_sample < N_t
                        image_tdma_2_transmit(yy,xx) = image_tdma_2_transmit(yy,xx) + m(t_sample);
                        if transmit == 1
                            image_tdma_1_transmit(yy,xx) = image_tdma_1_transmit(yy,xx) + m(t_sample);
                        end % if
                    end % if
                end % yy
            end % xx
        end % receiver
    end % transmit
    save('image_tdma_1_transmit.mat', 'image_tdma_1_transmit');
    save('image_tdma_2_transmit.mat', 'image_tdma_2_transmit');
end % if/else

%% TDMA two transmit
ploting(image_tdma_2_transmit, x, y, 'TDMA DAS, 2 transmits', dr)
print('images/tdma_2t',filetype)

%% TDMA one transmit
% dB
ploting(image_tdma_1_transmit, x, y, 'TDMA DAS, 1 transmit', dr)
print('images/tdma_1t',filetype)
%% Position of reflector
[max_value, max_index] = max((abs(image_tdma_2_transmit)),[],'all', 'linear');
[row,col] = ind2sub(size(image_tdma_2_transmit), max_index);
mpp = 1 / (xres / 10); % meter per pixel

% Cartesian
x_meter = (col * mpp) - 5
y_meter = row * mpp
% Polar
r = sqrt( (x_meter^2) + (y_meter^2) )
theta = atan2d(y_meter, x_meter)

%% Mesaured res.
% Lateral
k = 2000;
angle_range = 10;
angles = linspace((theta - angle_range/2),(theta + angle_range/2), k);
lateral_values = complex(zeros(1, k));
for i = [1:1:k]
    x_pos = (r * (cosd(angles(i)))); % meter
    y_pos = (r * (sind(angles(i)))); % meter
    row = round(y_pos / mpp); % row number
    col = round((x_pos + 5) / mpp); % column number
    lateral_values(i) = image_tdma_2_transmit(row, col);
end
max_value = max((abs(lateral_values)),[],'all');
normalized_lateral_values = (abs(lateral_values)) ./ max_value;
%figure()
%plot(angles, db(abs(normalized_lateral_values)));

tmp = normalized_lateral_values((db(abs(normalized_lateral_values))) >= -3); % FWHM part of curve
lateral_FWHM_deg = length(tmp) * (angle_range/k); % FWHM amount of degres
lateral_FWHM = deg2rad(lateral_FWHM_deg) * r


% Axial
k = 2000;
range_range = 0.2; % meter(s)
ranges = linspace((r - range_range/2),(r + range_range/2), k);
axial_values = complex(zeros(1, k));
for i = [1:1:k]
    x_pos = (ranges(i) * (cosd(theta))); % meter
    y_pos = (ranges(i) * (sind(theta))); % meter
    row = round(y_pos / mpp); % row number
    col = round((x_pos + 5) / mpp); % column number
    axial_values(i) = image_tdma_2_transmit(row, col);
end
max_value = max((abs(axial_values)),[],'all');
normalized_axial_values = (abs(axial_values)) ./ max_value;
%figure()
%plot(ranges, db(abs(normalized_axial_values)));
tmp = normalized_axial_values((db(abs(normalized_axial_values))) >= -3); % FWHM part of curve
axial_FWHM = length(tmp) * (range_range/k)

%% DAS CDMA
% xcorr over all data
all_m = zeros(N_t, N_tx, N_rx);
all_tm = zeros(N_t, N_tx, N_rx);

for receiver = 1:N_rx
    tmp_data = cdma_data(:,receiver); % single receve, both transmit
    [m1,lag1] = xcorr(tmp_data, LFM_down); % left transmiter 
    [m2,lag2] = xcorr(tmp_data, LFM_up); % right transmiter 
    m1 = m1(lag1 >= 0);
    m2 = m2(lag2 >= 0);
    %tm1 = lag1(lag1 >= 0)/fs;
    %tm2 = lag2(lag2 >= 0)/fs;
    all_m(:,1,receiver) = m1;
    all_m(:,2,receiver) = m2;
    %all_tm(:,1,receiver) = tm;
    %all_tm(:,2,receiver) = tm;
end


% DAS
x = linspace(-5,5,xres);
y = linspace(0,5,yres);
image_cdma = complex(zeros(yres,xres));

%load_data = 1; % 0 = load, 1 = calculate

if load_data == 0
    load('image_cdma.mat');
elseif load_data == 1
    for transmit = 1:N_tx
        for receiver = 1:N_rx
            m = all_m(:,transmit,receiver);
            for xx = 1:xres
                for yy = 1:yres
                    r_t = sqrt( (tx_pos(transmit)-x(xx))^2 + (y(yy))^2 ); 
                    r_r = sqrt( (rx_pos(receiver)-x(xx))^2 + (y(yy))^2 );
                    r = r_t + r_r;
                    t_delay = r / c;
                    t_sample = round(t_delay * fs);
                    if t_sample > 0 && t_sample < N_t
                        image_cdma(yy,xx) = image_cdma(yy,xx) + m(t_sample);    
                    end % if
                end % yy
            end % xx
        end % receiver
    end % transmit    
    save('image_cdma.mat', 'image_cdma');
end % if/elseif
%% CDMA Plot
% DB
ploting(image_cdma, x, y, 'CDMA DAS', dr)
print('images/cdma',filetype)

%% Calculations

% Approx lateral resolution at 4 m: L = R*a = 4 (m) * 0.0635 = 0.254 m
% Axial resolution = c / (2*B) = 0.0170 m

% X is -5 to 5 meters so 10 meters 
% Y is 0 to 5 so 5 meters
% 2*yres = xres
% axial pixels per meter = lateral pixlels per meter

% yres [pixels]/ 5 [m] -- xres [pixels] / 10 [m] 
% 1000 / 5 -- 2000 / 10 = 200 [pixels/meter]

% So 1 pixel is 1/200 meters -> 0.005 meters per pixel.
% 0.017 / 0.005 = 3.4 times the axial resoulution.
% 0.1012 / 0.005 = 20.2400 times latteral resolution at point.


%% TDMA with tapering

% xcorr over all data
all_m = zeros(N_t, N_tx, N_rx);
%all_tm = zeros(N_t, N_tx, N_rx);

LFM_hamming = hamming(length(LFM_up)).';
LFM_up_tapered = LFM_hamming .* LFM_up;

for transmit = 1:N_tx
    for receiver = 1:N_rx
        tmp_data = tdma_data(:,receiver,transmit); % single receve, single transmit
        receve_hamming = hamming(length(tmp_data));
        tmp_data_tapered = receve_hamming .* tmp_data;
        [m,lag] = xcorr(tmp_data_tapered, LFM_up_tapered);
        m = m(lag >= 0);
        tm = lag(lag >= 0)/fs;
        all_m(:,transmit,receiver) = m;
        %all_tm(:,transmit,receiver) = tm;
    end
end

% DAS
x = linspace(-5,5,xres);
y = linspace(0,5,yres);
image_tdma_2_transmit_tapered = complex(zeros(yres,xres));
image_tdma_1_transmit_tapered = complex(zeros(yres,xres));

%load_data = 1; % 0 = load, 1 = calculate

if load_data == 0
    load('image_tdma_1_transmit_tapered.mat');
    load('image_tdma_2_transmit_tapered.mat');
elseif load_data == 1
    for transmit = 1:N_tx
        for receiver = 1:N_rx
            m = all_m(:,transmit,receiver);
            for xx = 1:xres
                for yy = 1:yres
                    r_t = sqrt( (tx_pos(transmit)-x(xx))^2 + (y(yy))^2 ); 
                    r_r = sqrt( (rx_pos(receiver)-x(xx))^2 + (y(yy))^2 );
                    r = r_t + r_r;
                    t_delay = r / c;
                    t_sample = round(t_delay * fs);
                    if t_sample > 0 && t_sample < N_t
                        image_tdma_2_transmit_tapered(yy,xx) = image_tdma_2_transmit_tapered(yy,xx) + m(t_sample);
                        if transmit == 1
                            image_tdma_1_transmit_tapered(yy,xx) = image_tdma_1_transmit_tapered(yy,xx) + m(t_sample);
                        end % if
                    end % if
                end % yy
            end % xx
        end % receiver
    end % transmit
    save('image_tdma_1_transmit_tapered.mat', 'image_tdma_1_transmit_tapered');
    save('image_tdma_2_transmit_tapered.mat', 'image_tdma_2_transmit_tapered');
end % if/else

%% TDMA two transmit
ploting(image_tdma_2_transmit_tapered, x, y, 'TDMA tapered, 2 transmits', dr)
print('images/tdma_2t_tapered',filetype)

%% TDMA one transmit
% dB
ploting(image_tdma_1_transmit_tapered, x, y, 'TDMA tapered, 1 transmits', dr)
print('images/tdma_1t_tapered',filetype)

%% TDMA with virtual array

% xcorr over all data
all_m = zeros(N_t, N_tx, N_rx);
%all_tm = zeros(N_t, N_tx, N_rx);

for transmit = 1:N_tx
    for receiver = 1:N_rx
        tmp_data = tdma_data(:,receiver,transmit); % single receve, single transmit
        [m,lag] = xcorr(tmp_data, LFM_up);
        m = m(lag >= 0);
        tm = lag(lag >= 0)/fs;
        all_m(:,transmit,receiver) = m;
        %all_tm(:,transmit,receiver) = tm;
    end
end

% DAS
x = linspace(-5,5,xres);
y = linspace(0,5,yres);
image_tdma_2_transmit_virtual = complex(zeros(yres,xres));
image_tdma_1_transmit_virtual = complex(zeros(yres,xres));

%load_data = 1; % 0 = load, 1 = calculate

if load_data == 0
    load('image_tdma_1_transmit_virtual.mat');
    load('image_tdma_2_transmit_virtual.mat');
elseif load_data == 1
    for transmit = 1:N_tx
        for receiver = 1:N_rx
            m = all_m(:,transmit,receiver);
            for xx = 1:xres
                for yy = 1:yres
                    vx_pos = (tx_pos(transmit) + rx_pos(receiver)) / 2; % virtual element posisjon
                    r = 2 * (sqrt( (vx_pos - x(xx))^2 + (y(yy))^2 ));
                    t_delay = r / c;
                    t_sample = round(t_delay * fs);
                    if t_sample > 0 && t_sample < N_t
                        image_tdma_2_transmit_virtual(yy,xx) = image_tdma_2_transmit_virtual(yy,xx) + m(t_sample);
                        if transmit == 1
                            image_tdma_1_transmit_virtual(yy,xx) = image_tdma_1_transmit_virtual(yy,xx) + m(t_sample);
                        end % if
                    end % if
                end % yy
            end % xx
        end % receiver
    end % transmit
    save('image_tdma_1_transmit_virtual.mat', 'image_tdma_1_transmit_virtual');
    save('image_tdma_2_transmit_virtual.mat', 'image_tdma_2_transmit_virtual');
end % if/elseif

%% TDMA viritual two transmit
ploting(image_tdma_2_transmit_virtual, x, y, 'TDMA virtual, 2 transmits', dr)
print('images/tdma_2t_virtual',filetype)

%% TDMA viriutal one transmit
% dB
ploting(image_tdma_1_transmit_virtual, x, y, 'TDMA virtual, 1 transmits', dr)
print('images/tdma_1t_virtual',filetype)





