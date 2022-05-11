clear all
run('generate_data.m')
filetype = '-depsc'; % eps with colour

%% Task 1
R = x*x' / N;
R_max = max(abs(R),[],'all');
figure()
imagesc(abs(R)/R_max)
a = colorbar;
a.Label.String = 'Inter-element correlation';
xlabel('Elment number')
ylabel('Elment number')
title('Normalized spacial correlation matrix of R')
print('images/R',filetype)

%% Task 2
d = 0.5;
kd = 2 * pi * d;
DOA = (-40:0.25:50);
L = length(DOA);
P_DAS = zeros(1, L);
ula_spacing = [0:(M-1)]';
for n = (1:L)
    phi = -kd*sin(DOA(n)*pi / 180);
    a = exp(1j*phi).^ula_spacing;
    P_DAS(n) = (a'*R*a) / M;
end
ploting(DOA,P_DAS,"Clasical spectrum",'power');
print('images/Classical',filetype)

%% Task 3
P_capton = zeros(1, L);
Rinv = inv(R);
for n = (1:L)
    phi = -kd*sin(DOA(n)*pi / 180); 
    a = exp(1j*phi).^ula_spacing;
    P_capton(n) = 1 / (a'*Rinv*a);
end
ploting(DOA,P_capton,"Capon's beamformer",'power'); % Minimum variance?
print('images/Capon',filetype)

%% Task 4
[V, D] = eig(R); % matrix of column eigenvectors, diagonal eigenvalue matrix
[dd, I] = sort(diag(D)); % Sorted eigenvalues as well as coresponding index values
dd = flipud(dd); % eigenvalues in desending order
V = V(:,flipud(I)); % sorted columns by decreasing eigenvalue

figure()
plot(dd,'*','MarkerSize',15)
xlim('padded')
ylim('padded')
xlabel('Eigenvalue #')
ylabel('Eigenvalue')
title('Eigenvalue distribution')
print('images/Eigenvalues',filetype)

%% Task 5
sourses = 2;
P_MUSIC = zeros(1, L);
U = V(:,sourses+1:end); % Remove the signal+nose columns, leaving only noise coulumns 
PI = U*U'; % Orthogonal projector onto the noise subspace
for n = (1:L)
    phi = -kd*sin(DOA(n)*pi / 180); 
    a = exp(1j*phi).^ula_spacing;
    P_MUSIC(n) = (a'*a) / (a'*PI*a);
end
ploting(DOA,P_MUSIC,"Music beamformer",'power');
print('images/Music',filetype)

%% Task 6
sourses = 2;
P_EV = zeros(1, L);
U = V(:,sourses+1:end); % Remove the signal+nose columns, leaving only noise coulumns
A = inv(diag(dd(sourses+1:end)));
UAU = U*A*U';
for n = (1:L)
    phi = -kd*sin(DOA(n)*pi / 180);
    a = exp(1j*phi).^ula_spacing;
    P_EV(n) = (a'*a) / (a'*UAU*a);
end
ploting(DOA,P_MUSIC,"Eigenvector beamformer",'power');
print('images/EV',filetype)

%% Task 7
for sourses = [0,1,3]
    P_MUSIC = zeros(1, L);
    P_EV = zeros(1, L);
    U = V(:,sourses+1:end); % Remove the signal+nose columns, leaving only noise coulumns 
    PI = U*U'; % Orthogonal projector onto the noise subspace
    A = inv(diag(dd(sourses+1:end)));
    
    for n = (1:L)
        phi = -kd*sin(DOA(n)*pi / 180); 
        a = exp(1j*phi).^ula_spacing;
        P_MUSIC(n) = (a'*a) / (a'*PI*a);
        P_EV(n) = (a'*a) / (a'*U*A*U'*a);
    end
    
    format1 = "Music beamformer (sources = %i)";
    format2 = 'images/Music_%i_sources';
    ploting(DOA,P_MUSIC,sprintf(format1,sourses),'power');
    if sourses == 0
        ylim([-0.001 0.001])
    end
    print(sprintf(format2, sourses),filetype)
    
    format3 = "Eigenvector beamformer (sources = %i)";
    format4 = 'images/EV_%i_sources';
    ploting(DOA,P_EV,sprintf(format3,sourses),'power');
    print(sprintf(format4,sourses),filetype)
end

