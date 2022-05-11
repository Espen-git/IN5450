clear all
diary on

P = ASAparameters;
lmd = P.lambda;
d = lmd/2;
L = 1001;
NEls = 101;
NPos = 25;
weights = ones(1,NPos)*(1/NPos);
theta = linspace(-pi/2,pi/2,L);
u = sin(theta);
kx = 2*pi/lmd*(u);

anas = [];
figure();
for i = 1:100
    pos = [];
    while length(pos)<NPos - 2
    pos = unique(ceil((NEls-2)*rand(1,NPos*2)),'stable');
    end
    ElPos = [-(NEls-1)/2 sort(pos(1:NPos-2))-(NEls-1)/2 (NEls-1)/2]*d;

    Resp = beampattern(ElPos, kx, weights);
    ana = analyzeBP(u, Resp);
    anas = [anas, ana];

    subplot(2,1,1)
    plot1 = plot(u,db(abs(Resp)/max(abs(Resp))),'k'); grid on; ylim([-30 0]);
    plot1.Color(4) = 0.075;
    title(['Plot of several thinned arrays (', ...
    int2str(NPos),' of ',int2str(NEls),' elements)']);
    ylabel('dB');
    xlabel('u=sin $\theta$','interpreter','latex');
    hold on
    subplot(2,1,2)
    plot2 = plot([ElPos; ElPos],[zeros(1,NPos); ones(1,NPos)],'k');
    for ii=1:length(plot2)
    plot2(ii).Color(4) = 0.075;
    end
    title(['Distribution of thinned elements']);
    ylabel('El weight');
    xlabel('Elementposition');
    hold on
end

values = struct2cell(anas);
values = cell2mat(values);
ML_width_values = reshape(values(4,:,:),[1,100]);
SL_max_values = reshape(values(3,:,:),[1,100]);

ML_width_mean = mean(ML_width_values)
ML_width_sd = std(ML_width_values)
ML_width_max = max(ML_width_values);

SL_mean = mean(SL_max_values)
SL_sd = std(SL_max_values);
SL_max = max(SL_max_values)

%%
% 'Optimal' arrays
elements_binary = [1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 1 1 0 0 0 1 0 0 1 0 0 1 0 1 0 1 1 1 1 1 0 0 1 0 0 0 1 0 0 0 0 0 1 0 0 1 0 0 1 0 0 0 0 0 0 0 1 1 0 0 0 1]; 
ElPos = find(elements_binary)*d;

Resp = beampattern(ElPos, kx, weights);
ana = analyzeBP(u, Resp)

figure();
plot(u, db(abs(Resp)))
ylim([-30 0]);
xlabel('u [sin(theta)]'); ylabel('[dB]');

%%
% N=50

NPos = 50;
weights = ones(1,NPos)*(1/NPos);

anas = [];
figure();
for i = 1:100
    pos = [];
    while length(pos)<NPos - 2
    pos = unique(ceil((NEls-2)*rand(1,NPos*2)),'stable');
    end
    ElPos = [-(NEls-1)/2 sort(pos(1:NPos-2))-(NEls-1)/2 (NEls-1)/2]*d;

    Resp = beampattern(ElPos, kx, weights);
    ana = analyzeBP(u, Resp);
    anas = [anas, ana];

    subplot(2,1,1)
    plot1 = plot(u,db(abs(Resp)/max(abs(Resp))),'k'); grid on; ylim([-30 0]);
    plot1.Color(4) = 0.075;
    title(['Plot of several thinned arrays (', ...
    int2str(NPos),' of ',int2str(NEls),' elements)']);
    ylabel('dB');
    xlabel('u=sin $\theta$','interpreter','latex');
    hold on
    subplot(2,1,2)
    plot2 = plot([ElPos; ElPos],[zeros(1,NPos); ones(1,NPos)],'k');
    for ii=1:length(plot2)
    plot2(ii).Color(4) = 0.075;
    end
    title(['Distribution of thinned elements']);
    ylabel('El weight');
    xlabel('Elementposition');
    hold on
end

values = struct2cell(anas);
values = cell2mat(values);
ML_width_values = reshape(values(4,:,:),[1,100]);
SL_max_values = reshape(values(3,:,:),[1,100]);

ML_width_mean = mean(ML_width_values)
ML_width_sd = std(ML_width_values)
ML_width_max = max(ML_width_values);

SL_mean = mean(SL_max_values)
SL_sd = std(SL_max_values);
SL_max = max(SL_max_values)

%%
% N=75

NPos = 75;
weights = ones(1,NPos)*(1/NPos);
theta = linspace(-pi/2,pi/2,L);
u = sin(theta);
kx = 2*pi/lmd*(u);

anas = [];
figure();
for i = 1:100
    pos = [];
    while length(pos)<NPos - 2
    pos = unique(ceil((NEls-2)*rand(1,NPos*2)),'stable');
    end
    ElPos = [-(NEls-1)/2 sort(pos(1:NPos-2))-(NEls-1)/2 (NEls-1)/2]*d;

    Resp = beampattern(ElPos, kx, weights);
    ana = analyzeBP(u, Resp);
    anas = [anas, ana];

    subplot(2,1,1)
    plot1 = plot(u,db(abs(Resp)/max(abs(Resp))),'k'); grid on; ylim([-30 0]);
    plot1.Color(4) = 0.075;
    title(['Plot of several thinned arrays (', ...
    int2str(NPos),' of ',int2str(NEls),' elements)']);
    ylabel('dB');
    xlabel('u=sin $\theta$','interpreter','latex');
    hold on
    subplot(2,1,2)
    plot2 = plot([ElPos; ElPos],[zeros(1,NPos); ones(1,NPos)],'k');
    for ii=1:length(plot2)
    plot2(ii).Color(4) = 0.075;
    end
    title(['Distribution of thinned elements']);
    ylabel('El weight');
    xlabel('Elementposition');
    hold on
end

values = struct2cell(anas);
values = cell2mat(values);
ML_width_values = reshape(values(4,:,:),[1,100]);
SL_max_values = reshape(values(3,:,:),[1,100]);

ML_width_mean = mean(ML_width_values)
ML_width_sd = std(ML_width_values)
ML_width_max = max(ML_width_values);

SL_mean = mean(SL_max_values)
SL_sd = std(SL_max_values);
SL_max = max(SL_max_values)













