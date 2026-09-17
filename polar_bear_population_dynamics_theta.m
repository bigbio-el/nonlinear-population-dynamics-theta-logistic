% Polar Bears Population Dynamics Analysis
% Authors: Evgeniia Lavrenteva, Michael Binns
% Revised Chukchi Sea / Alaska-Chukotka population model

clear;
close all;
clc;


%% Model parameters

% Chukchi Sea population size estimate
N0 = 2937;
N0_low = 1552;
N0_high = 5944;

% Lower-growth demographic scenario
rMNPL1 = 0.02;
rMNPL1_low = 0.00;
rMNPL1_high = 0.06;

xMNPL1 = 0.73;
xMNPL1_low = 0.69;
xMNPL1_high = 0.74;

% Higher-growth demographic scenario
rMNPL2 = 0.04;
rMNPL2_low = 0.00;
rMNPL2_high = 0.07;

xMNPL2 = 0.70;
xMNPL2_low = 0.68;
xMNPL2_high = 0.74;

% Harvesting scenarios, bears per year
H_recent = 44.4;
H_historical = 56;
H_low = 50;
H_management = 85;
H_high = 120;

% Decline in environmental carrying capacity per decade
dK0 = 0.00;
dK1 = 0.03;
dK2 = 0.09;

% Numerical settings
dt = 0.05;
tspan = 0:dt:35;


%% Theta-logistic parameters

% Calculate theta from relative density

fun1 = @(theta) (1/(theta+1))^(1/theta) - xMNPL1;
theta1 = fzero(fun1,[0.01 50]);

fun2 = @(theta) (1/(theta+1))^(1/theta) - xMNPL2;
theta2 = fzero(fun2,[0.01 50]);

% Maximum intrinsic per-capita growth rate
r1 = rMNPL1/(1-xMNPL1^theta1);
r2 = rMNPL2/(1-xMNPL2^theta2);

% Reference carrying-capacity scale
K1 = N0/xMNPL1;
K2 = N0/xMNPL2;

fprintf('\nCalculated model parameters\n');
fprintf('Lower-growth: theta = %.4f, rmax = %.4f yr^-1, K0 = %.1f bears\n', ...
    theta1,r1,K1);
fprintf('Higher-growth: theta = %.4f, rmax = %.4f yr^-1, K0 = %.1f bears\n', ...
    theta2,r2,K2);
fprintf('rMNPL/rmax = %.3f and %.3f\n\n',rMNPL1/r1,rMNPL2/r2);


%% Production maximum

Pmax1_0 = rMNPL1*N0;
Pmax2_0 = rMNPL2*N0;

Pmax1_K1_35 = Pmax1_0*(1-dK1)^(35/10);
Pmax1_K2_35 = Pmax1_0*(1-dK2)^(35/10);
Pmax2_K1_35 = Pmax2_0*(1-dK1)^(35/10);
Pmax2_K2_35 = Pmax2_0*(1-dK2)^(35/10);

% Crossing time of the production maximum and reference harvest
tcross1_K1 = 10*log(H_recent/Pmax1_0)/log(1-dK1);
tcross1_K2 = 10*log(H_recent/Pmax1_0)/log(1-dK2);
tcross2_K1 = 10*log(H_recent/Pmax2_0)/log(1-dK1);
tcross2_K2 = 10*log(H_recent/Pmax2_0)/log(1-dK2);


%% Plot settings

blue = [0 114 178]/255;
orange = [213 94 0]/255;
green = [0 158 115]/255;
sky = [86 180 233]/255;
amber = [230 159 0]/255;
charcoal = [0.18 0.18 0.18];
midgray = [0.45 0.45 0.45];
lightgray = [0.86 0.86 0.86];
verylight = [0.93 0.93 0.93];

font_name = 'Helvetica';
font_size = 9;
label_size = 10;
title_size = 10;
legend_size = 8;
line_width = 1.8;
line_width2 = 2.2;
marker_size = 6;

outdir = 'figures';
if ~exist(outdir,'dir')
    mkdir(outdir);
end


%% Figure 2

density = 0:0.002:1.2;

production1 = r1*density.*(1-density.^theta1);
production2 = r2*density.*(1-density.^theta2);

% Published Chukchi Sea population growth estimate
g_observed = -0.01;
g_observed_low = -0.04;
g_observed_high = 0.03;

% Initial net growth under the historical harvesting estimate
g_model1 = rMNPL1-H_historical/N0;
g_model2 = rMNPL2-H_historical/N0;

f = figure('Color','w','Units','centimeters','Position',[2 2 18 8.2]);
set(f,'Renderer','painters');

subplot(1,2,1);
plot(density,production1,'Color',blue,'LineStyle','-','LineWidth',line_width2); hold on;
plot(density,production2,'Color',orange,'LineStyle','--','LineWidth',line_width2);
plot(xMNPL1,xMNPL1*rMNPL1,'o','Color',blue,'MarkerFaceColor',blue, ...
    'MarkerSize',marker_size,'LineWidth',1.0);
plot(xMNPL2,xMNPL2*rMNPL2,'s','Color',orange,'MarkerFaceColor',orange, ...
    'MarkerSize',marker_size,'LineWidth',1.0);
plot([0 1.2],[0 0],'-','Color',lightgray,'LineWidth',0.9);
hold off;
xlabel('Relative density, N/K','FontSize',label_size);
ylabel('Relative net population production','FontSize',label_size);
xlim([0 1.2]); ylim([-0.10 0.04]);
legend('Lower-growth','Higher-growth','MNPL, lower','MNPL, higher', ...
    'Location','southwest','FontSize',legend_size,'Box','off');
text(0.02,0.96,'(a)','Units','normalized','FontWeight','bold','FontSize',title_size, ...
    'VerticalAlignment','top');

subplot(1,2,2);
errorbar(1,g_observed,g_observed-g_observed_low,g_observed_high-g_observed, ...
    'o','Color',charcoal,'MarkerFaceColor','w','MarkerSize',marker_size+1, ...
    'LineWidth',1.3,'CapSize',8); hold on;
plot(2,g_model1,'s','Color',blue,'MarkerFaceColor',blue,'MarkerSize',marker_size+1);
plot(3,g_model2,'^','Color',orange,'MarkerFaceColor',orange,'MarkerSize',marker_size+1);
plot([0.55 3.45],[0 0],'--','Color',midgray,'LineWidth',1.0);
hold off;
ylabel('Net population growth rate (yr^{-1})','FontSize',label_size);
set(gca,'XTick',[1 2 3],'XTickLabel',{'Published','Lower-growth','Higher-growth'});
xlim([0.55 3.45]); ylim([-0.045 0.035]);
text(0.02,0.96,'(b)','Units','normalized','FontWeight','bold','FontSize',title_size, ...
    'VerticalAlignment','top');

axs = findall(f,'Type','axes');
set(axs,'FontName',font_name,'FontSize',font_size,'LineWidth',0.8, ...
    'TickDir','out','Box','off','XColor',charcoal,'YColor',charcoal);
for k = 1:length(axs)
    grid(axs(k),'on');
    axs(k).GridColor = lightgray;
    axs(k).GridAlpha = 0.45;
end

drawnow;
savefig(f,fullfile(outdir,'Figure2_density_consistency.fig'));
set(f,'PaperPositionMode','auto');
print(f,fullfile(outdir,'Figure2_density_consistency.png'),'-dpng','-r600');


%% Figure 3
Hvalues = [0 H_recent H_low H_management H_high];
dKvalues = [dK0 dK1 dK2];

N_harvest1 = zeros(length(Hvalues),length(tspan));
N_harvest2 = zeros(length(Hvalues),length(tspan));

% Lower-growth scenario
for j = 1:length(Hvalues)

    H = Hvalues(j);
    N = zeros(size(tspan));
    N(1) = N0;

    for i = 2:length(tspan)

        t = tspan(i-1);

        odefun = @(tt,y) r1*max(y,0)*(1-(max(y,0)/(K1*(1-dK1)^(tt/10)))^theta1)-H;

        k1 = odefun(t,N(i-1));
        k2 = odefun(t+dt/2,N(i-1)+dt*k1/2);
        k3 = odefun(t+dt/2,N(i-1)+dt*k2/2);
        k4 = odefun(t+dt,N(i-1)+dt*k3);

        y = N(i-1)+dt*(k1+2*k2+2*k3+k4)/6;

        if y < 0
            y = 0;
        end

        N(i) = y;

    end

    N_harvest1(j,:) = N;

end

% Higher-growth scenario
for j = 1:length(Hvalues)

    H = Hvalues(j);
    N = zeros(size(tspan));
    N(1) = N0;

    for i = 2:length(tspan)

        t = tspan(i-1);

        odefun = @(tt,y) r2*max(y,0)*(1-(max(y,0)/(K2*(1-dK1)^(tt/10)))^theta2)-H;

        k1 = odefun(t,N(i-1));
        k2 = odefun(t+dt/2,N(i-1)+dt*k1/2);
        k3 = odefun(t+dt/2,N(i-1)+dt*k2/2);
        k4 = odefun(t+dt,N(i-1)+dt*k3);

        y = N(i-1)+dt*(k1+2*k2+2*k3+k4)/6;

        if y < 0
            y = 0;
        end

        N(i) = y;

    end

    N_harvest2(j,:) = N;

end

N_habitat1 = zeros(length(dKvalues),length(tspan));
N_habitat2 = zeros(length(dKvalues),length(tspan));

% Lower-growth scenario
for j = 1:length(dKvalues)

    dK = dKvalues(j);
    H = H_recent;
    N = zeros(size(tspan));
    N(1) = N0;

    for i = 2:length(tspan)

        t = tspan(i-1);

        odefun = @(tt,y) r1*max(y,0)*(1-(max(y,0)/(K1*(1-dK)^(tt/10)))^theta1)-H;

        k1 = odefun(t,N(i-1));
        k2 = odefun(t+dt/2,N(i-1)+dt*k1/2);
        k3 = odefun(t+dt/2,N(i-1)+dt*k2/2);
        k4 = odefun(t+dt,N(i-1)+dt*k3);

        y = N(i-1)+dt*(k1+2*k2+2*k3+k4)/6;

        if y < 0
            y = 0;
        end

        N(i) = y;

    end

    N_habitat1(j,:) = N;

end

% Higher-growth scenario
for j = 1:length(dKvalues)

    dK = dKvalues(j);
    H = H_recent;
    N = zeros(size(tspan));
    N(1) = N0;

    for i = 2:length(tspan)

        t = tspan(i-1);

        odefun = @(tt,y) r2*max(y,0)*(1-(max(y,0)/(K2*(1-dK)^(tt/10)))^theta2)-H;

        k1 = odefun(t,N(i-1));
        k2 = odefun(t+dt/2,N(i-1)+dt*k1/2);
        k3 = odefun(t+dt/2,N(i-1)+dt*k2/2);
        k4 = odefun(t+dt,N(i-1)+dt*k3);

        y = N(i-1)+dt*(k1+2*k2+2*k3+k4)/6;

        if y < 0
            y = 0;
        end

        N(i) = y;

    end

    N_habitat2(j,:) = N;

end

f = figure('Color','w','Units','centimeters','Position',[2 2 18 13.2]);
set(f,'Renderer','painters');

harvestColors = [charcoal; green; sky; amber; orange];
harvestStyles = {'-','--',':','-.','-'};
harvestWidths = [1.8 2.2 2.0 2.0 2.1];

harvestLabels = {sprintf('%.0f',Hvalues(1)),sprintf('%.1f',Hvalues(2)), ...
    sprintf('%.0f',Hvalues(3)),sprintf('%.0f',Hvalues(4)),sprintf('%.0f',Hvalues(5))};
habitatLabels = {'Stable K',sprintf('%.0f%% per decade',100*dK1), ...
    sprintf('%.0f%% per decade',100*dK2)};

subplot(2,2,1); hold on;
for j = 1:length(Hvalues)
    plot(tspan,N_harvest1(j,:),'Color',harvestColors(j,:), ...
        'LineStyle',harvestStyles{j},'LineWidth',harvestWidths(j));
end
plot([tspan(1) tspan(end)],[N0 N0],'--','Color',midgray,'LineWidth',1.0);
hold off;
ylabel('Population size','FontSize',label_size);
xlim([tspan(1) tspan(end)]); ylim([0 4100]);
legend(harvestLabels,'Location','southwest', ...
    'FontSize',legend_size,'Box','off');
title('Lower-growth','FontSize',title_size,'FontWeight','bold');
text(0.02,0.96,'(a)','Units','normalized','FontWeight','bold','FontSize',title_size, ...
    'VerticalAlignment','top');

subplot(2,2,2); hold on;
for j = 1:length(Hvalues)
    plot(tspan,N_harvest2(j,:),'Color',harvestColors(j,:), ...
        'LineStyle',harvestStyles{j},'LineWidth',harvestWidths(j));
end
plot([tspan(1) tspan(end)],[N0 N0],'--','Color',midgray,'LineWidth',1.0);
hold off;
xlim([tspan(1) tspan(end)]); ylim([0 4100]);
legend(harvestLabels,'Location','southwest', ...
    'FontSize',legend_size,'Box','off');
title('Higher-growth','FontSize',title_size,'FontWeight','bold');
text(0.02,0.96,'(b)','Units','normalized','FontWeight','bold','FontSize',title_size, ...
    'VerticalAlignment','top');

subplot(2,2,3); hold on;
plot(tspan,N_habitat1(1,:),'Color',charcoal,'LineStyle','-','LineWidth',line_width);
plot(tspan,N_habitat1(2,:),'Color',green,'LineStyle','--','LineWidth',line_width2);
plot(tspan,N_habitat1(3,:),'Color',orange,'LineStyle',':','LineWidth',line_width2);
plot([tspan(1) tspan(end)],[N0 N0],'--','Color',midgray,'LineWidth',1.0);
hold off;
xlabel('Years','FontSize',label_size);
ylabel('Population size','FontSize',label_size);
xlim([tspan(1) tspan(end)]); ylim([2400 4100]);
legend(habitatLabels,'Location','northwest', ...
    'FontSize',legend_size,'Box','off');
title('Lower-growth','FontSize',title_size,'FontWeight','bold');
text(0.02,0.96,'(c)','Units','normalized','FontWeight','bold','FontSize',title_size, ...
    'VerticalAlignment','top');

subplot(2,2,4); hold on;
plot(tspan,N_habitat2(1,:),'Color',charcoal,'LineStyle','-','LineWidth',line_width);
plot(tspan,N_habitat2(2,:),'Color',green,'LineStyle','--','LineWidth',line_width2);
plot(tspan,N_habitat2(3,:),'Color',orange,'LineStyle',':','LineWidth',line_width2);
plot([tspan(1) tspan(end)],[N0 N0],'--','Color',midgray,'LineWidth',1.0);
hold off;
xlabel('Years','FontSize',label_size);
xlim([tspan(1) tspan(end)]); ylim([2400 4100]);
legend(habitatLabels,'Location','southeast', ...
    'FontSize',legend_size,'Box','off');
title('Higher-growth','FontSize',title_size,'FontWeight','bold');
text(0.02,0.96,'(d)','Units','normalized','FontWeight','bold','FontSize',title_size, ...
    'VerticalAlignment','top');

axs = findall(f,'Type','axes');
set(axs,'FontName',font_name,'FontSize',font_size,'LineWidth',0.8, ...
    'TickDir','out','Box','off','XColor',charcoal,'YColor',charcoal);
for k = 1:length(axs)
    grid(axs(k),'on');
    axs(k).GridColor = lightgray;
    axs(k).GridAlpha = 0.42;
end

drawnow;
savefig(f,fullfile(outdir,'Figure3_population_trajectories.fig'));
set(f,'PaperPositionMode','auto');
print(f,fullfile(outdir,'Figure3_population_trajectories.png'),'-dpng','-r600');


%% Figure 4

harvest_values = 0:2:H_high;
decline_values = 0:0.002:dK2;

final_population1 = zeros(length(decline_values),length(harvest_values));
final_population2 = zeros(length(decline_values),length(harvest_values));

% Lower-growth scenario
for a = 1:length(decline_values)

    dK = decline_values(a);

    for b = 1:length(harvest_values)

        H = harvest_values(b);
        N = zeros(size(tspan));
        N(1) = N0;

        for i = 2:length(tspan)

            t = tspan(i-1);

            odefun = @(tt,y) r1*max(y,0)*(1-(max(y,0)/(K1*(1-dK)^(tt/10)))^theta1)-H;

            k1 = odefun(t,N(i-1));
            k2 = odefun(t+dt/2,N(i-1)+dt*k1/2);
            k3 = odefun(t+dt/2,N(i-1)+dt*k2/2);
            k4 = odefun(t+dt,N(i-1)+dt*k3);

            y = N(i-1)+dt*(k1+2*k2+2*k3+k4)/6;

            if y < 0
                y = 0;
            end

            N(i) = y;

        end

        final_population1(a,b) = N(end);

    end

end

% Higher-growth scenario
for a = 1:length(decline_values)

    dK = decline_values(a);

    for b = 1:length(harvest_values)

        H = harvest_values(b);
        N = zeros(size(tspan));
        N(1) = N0;

        for i = 2:length(tspan)

            t = tspan(i-1);

            odefun = @(tt,y) r2*max(y,0)*(1-(max(y,0)/(K2*(1-dK)^(tt/10)))^theta2)-H;

            k1 = odefun(t,N(i-1));
            k2 = odefun(t+dt/2,N(i-1)+dt*k1/2);
            k3 = odefun(t+dt/2,N(i-1)+dt*k2/2);
            k4 = odefun(t+dt,N(i-1)+dt*k3);

            y = N(i-1)+dt*(k1+2*k2+2*k3+k4)/6;

            if y < 0
                y = 0;
            end

            N(i) = y;

        end

        final_population2(a,b) = N(end);

    end

end

remaining1 = 100*final_population1/N0;
remaining2 = 100*final_population2/N0;

f = figure('Color','w','Units','centimeters','Position',[2 2 17.2 7.9]);
set(f,'Renderer','painters');

tl = tiledlayout(f,1,2,'TileSpacing','compact','Padding','compact');

contour_width = 1.25;
reference_marker_size = 5.5;
reference_marker_width = 0.9;

% (a) Lower-growth scenario
ax1 = nexttile(tl,1);
imagesc(ax1,harvest_values,decline_values*100,remaining1);
set(ax1,'YDir','normal');
hold(ax1,'on');

contour(ax1,harvest_values,decline_values*100,remaining1,[100 100], ...
    'Color',charcoal,'LineWidth',contour_width);

plot(ax1,H_recent,100*dK1,'o','MarkerEdgeColor',charcoal, ...
    'MarkerFaceColor','w','MarkerSize',reference_marker_size, ...
    'LineWidth',reference_marker_width);
plot(ax1,H_management,100*dK1,'s','MarkerEdgeColor',charcoal, ...
    'MarkerFaceColor','w','MarkerSize',reference_marker_size, ...
    'LineWidth',reference_marker_width);
plot(ax1,H_recent,100*dK2,'^','MarkerEdgeColor',charcoal, ...
    'MarkerFaceColor','w','MarkerSize',reference_marker_size, ...
    'LineWidth',reference_marker_width);

hold(ax1,'off');

title(ax1,'Lower-growth','FontName',font_name,'FontSize',title_size, ...
    'FontWeight','bold');
xlabel(ax1,'Annual harvest (bears yr^{-1})','FontName',font_name, ...
    'FontSize',label_size);
ylabel(ax1,'Decline in K (% per decade)','FontName',font_name, ...
    'FontSize',label_size);

text(ax1,0.025,0.965,'(a)','Units','normalized','FontName',font_name, ...
    'FontWeight','bold','FontSize',title_size,'VerticalAlignment','top', ...
    'Color',charcoal);

xlim(ax1,[0 H_high]);
ylim(ax1,[0 100*dK2+0.15]);
xticks(ax1,[0 30 60 90 120]);
yticks(ax1,0:1:9);
caxis(ax1,[0 150]);

% (b) Higher-growth scenario
ax2 = nexttile(tl,2);
imagesc(ax2,harvest_values,decline_values*100,remaining2);
set(ax2,'YDir','normal');
hold(ax2,'on');

contour(ax2,harvest_values,decline_values*100,remaining2,[100 100], ...
    'Color',charcoal,'LineWidth',contour_width);

plot(ax2,H_recent,100*dK1,'o','MarkerEdgeColor',charcoal, ...
    'MarkerFaceColor','w','MarkerSize',reference_marker_size, ...
    'LineWidth',reference_marker_width);
plot(ax2,H_management,100*dK1,'s','MarkerEdgeColor',charcoal, ...
    'MarkerFaceColor','w','MarkerSize',reference_marker_size, ...
    'LineWidth',reference_marker_width);
plot(ax2,H_recent,100*dK2,'^','MarkerEdgeColor',charcoal, ...
    'MarkerFaceColor','w','MarkerSize',reference_marker_size, ...
    'LineWidth',reference_marker_width);

hold(ax2,'off');

title(ax2,'Higher-growth','FontName',font_name,'FontSize',title_size, ...
    'FontWeight','bold');
xlabel(ax2,'Annual harvest (bears yr^{-1})','FontName',font_name, ...
    'FontSize',label_size);

text(ax2,0.025,0.965,'(b)','Units','normalized','FontName',font_name, ...
    'FontWeight','bold','FontSize',title_size,'VerticalAlignment','top', ...
    'Color',charcoal);

xlim(ax2,[0 H_high]);
ylim(ax2,[0 100*dK2+0.15]);
xticks(ax2,[0 30 60 90 120]);
yticks(ax2,0:1:9);
caxis(ax2,[0 150]);

colormap(f,parula(256));
set([ax1 ax2],'FontName',font_name,'FontSize',font_size,'LineWidth',0.8, ...
    'TickDir','out','Box','on','Layer','top', ...
    'XColor',charcoal,'YColor',charcoal);

ax1.GridColor = lightgray;
ax2.GridColor = lightgray;
ax1.GridAlpha = 0.18;
ax2.GridAlpha = 0.18;
ax1.XGrid = 'off'; ax1.YGrid = 'off';
ax2.XGrid = 'off'; ax2.YGrid = 'off';

cb = colorbar(ax2,'eastoutside');
cb.Layout.Tile = 'east';
cb.Label.String = 'N(35)/N_0 (%)';
cb.Label.FontName = font_name;
cb.Label.FontSize = label_size;
cb.FontName = font_name;
cb.FontSize = font_size;
cb.Ticks = [0 50 100 150];
cb.LineWidth = 0.8;
cb.Color = charcoal;

drawnow;
savefig(f,fullfile(outdir,'Figure4_habitat_harvest_surface.fig'));
set(f,'PaperPositionMode','auto');
print(f,fullfile(outdir,'Figure4_habitat_harvest_surface.png'),'-dpng','-r600');


%% Exact harvesting threshold, H*

Hstar1 = nan(size(decline_values));
Hstar2 = nan(size(decline_values));
threshold_tolerance = 0.01;

% Lower-growth scenario
for a = 1:length(decline_values)

    dK = decline_values(a);
    Hlo = 0;
    Hhi = H_high;

    H = Hlo;
    N = zeros(size(tspan));
    N(1) = N0;

    for i = 2:length(tspan)

        t = tspan(i-1);
        odefun = @(tt,y) r1*max(y,0)*(1-(max(y,0)/(K1*(1-dK)^(tt/10)))^theta1)-H;

        k1 = odefun(t,N(i-1));
        k2 = odefun(t+dt/2,N(i-1)+dt*k1/2);
        k3 = odefun(t+dt/2,N(i-1)+dt*k2/2);
        k4 = odefun(t+dt,N(i-1)+dt*k3);

        y = N(i-1)+dt*(k1+2*k2+2*k3+k4)/6;
        if y < 0
            y = 0;
        end
        N(i) = y;

    end

    Flo = N(end)-N0;

    H = Hhi;
    N = zeros(size(tspan));
    N(1) = N0;

    for i = 2:length(tspan)

        t = tspan(i-1);
        odefun = @(tt,y) r1*max(y,0)*(1-(max(y,0)/(K1*(1-dK)^(tt/10)))^theta1)-H;

        k1 = odefun(t,N(i-1));
        k2 = odefun(t+dt/2,N(i-1)+dt*k1/2);
        k3 = odefun(t+dt/2,N(i-1)+dt*k2/2);
        k4 = odefun(t+dt,N(i-1)+dt*k3);

        y = N(i-1)+dt*(k1+2*k2+2*k3+k4)/6;
        if y < 0
            y = 0;
        end
        N(i) = y;

    end

    Fhi = N(end)-N0;

    if Flo >= 0 && Fhi <= 0

        while Hhi-Hlo > threshold_tolerance

            Hmid = (Hlo+Hhi)/2;
            H = Hmid;
            N = zeros(size(tspan));
            N(1) = N0;

            for i = 2:length(tspan)

                t = tspan(i-1);
                odefun = @(tt,y) r1*max(y,0)*(1-(max(y,0)/(K1*(1-dK)^(tt/10)))^theta1)-H;

                k1 = odefun(t,N(i-1));
                k2 = odefun(t+dt/2,N(i-1)+dt*k1/2);
                k3 = odefun(t+dt/2,N(i-1)+dt*k2/2);
                k4 = odefun(t+dt,N(i-1)+dt*k3);

                y = N(i-1)+dt*(k1+2*k2+2*k3+k4)/6;
                if y < 0
                    y = 0;
                end
                N(i) = y;

            end

            Fmid = N(end)-N0;

            if Fmid >= 0
                Hlo = Hmid;
            else
                Hhi = Hmid;
            end

        end

        Hstar1(a) = Hlo;

    end

end

% Higher-growth scenario
for a = 1:length(decline_values)

    dK = decline_values(a);
    Hlo = 0;
    Hhi = H_high;

    H = Hlo;
    N = zeros(size(tspan));
    N(1) = N0;

    for i = 2:length(tspan)

        t = tspan(i-1);
        odefun = @(tt,y) r2*max(y,0)*(1-(max(y,0)/(K2*(1-dK)^(tt/10)))^theta2)-H;

        k1 = odefun(t,N(i-1));
        k2 = odefun(t+dt/2,N(i-1)+dt*k1/2);
        k3 = odefun(t+dt/2,N(i-1)+dt*k2/2);
        k4 = odefun(t+dt,N(i-1)+dt*k3);

        y = N(i-1)+dt*(k1+2*k2+2*k3+k4)/6;
        if y < 0
            y = 0;
        end
        N(i) = y;

    end

    Flo = N(end)-N0;

    H = Hhi;
    N = zeros(size(tspan));
    N(1) = N0;

    for i = 2:length(tspan)

        t = tspan(i-1);
        odefun = @(tt,y) r2*max(y,0)*(1-(max(y,0)/(K2*(1-dK)^(tt/10)))^theta2)-H;

        k1 = odefun(t,N(i-1));
        k2 = odefun(t+dt/2,N(i-1)+dt*k1/2);
        k3 = odefun(t+dt/2,N(i-1)+dt*k2/2);
        k4 = odefun(t+dt,N(i-1)+dt*k3);

        y = N(i-1)+dt*(k1+2*k2+2*k3+k4)/6;
        if y < 0
            y = 0;
        end
        N(i) = y;

    end

    Fhi = N(end)-N0;

    if Flo >= 0 && Fhi <= 0

        while Hhi-Hlo > threshold_tolerance

            Hmid = (Hlo+Hhi)/2;
            H = Hmid;
            N = zeros(size(tspan));
            N(1) = N0;

            for i = 2:length(tspan)

                t = tspan(i-1);
                odefun = @(tt,y) r2*max(y,0)*(1-(max(y,0)/(K2*(1-dK)^(tt/10)))^theta2)-H;

                k1 = odefun(t,N(i-1));
                k2 = odefun(t+dt/2,N(i-1)+dt*k1/2);
                k3 = odefun(t+dt/2,N(i-1)+dt*k2/2);
                k4 = odefun(t+dt,N(i-1)+dt*k3);

                y = N(i-1)+dt*(k1+2*k2+2*k3+k4)/6;
                if y < 0
                    y = 0;
                end
                N(i) = y;

            end

            Fmid = N(end)-N0;

            if Fmid >= 0
                Hlo = Hmid;
            else
                Hhi = Hmid;
            end

        end

        Hstar2(a) = Hlo;

    end

end


%% Figure 5

index_K1 = find(abs(decline_values-dK1)<1e-12,1);
index_K2 = find(abs(decline_values-dK2)<1e-12,1);

Hstar1_K1 = Hstar1(index_K1);
Hstar2_K1 = Hstar2(index_K1);
Hstar1_K2 = Hstar1(index_K2);
Hstar2_K2 = Hstar2(index_K2);

f = figure('Color','w','Units','centimeters','Position',[2 2 16.5 9.5]);
set(f,'Renderer','painters');

plot(decline_values*100,Hstar1,'Color',blue,'LineStyle','-','LineWidth',line_width2); hold on;
plot(decline_values*100,Hstar2,'Color',orange,'LineStyle','--','LineWidth',line_width2);
plot([0 100*dK2],[H_recent H_recent],'--','Color',green,'LineWidth',1.2);
plot([0 100*dK2],[H_management H_management],'-.','Color',amber,'LineWidth',1.2);
plot([100*dK1 100*dK1],[0 125],':','Color',midgray,'LineWidth',1.0);
plot([100*dK2 100*dK2],[0 125],':','Color',midgray,'LineWidth',1.0);
plot([100*dK1 100*dK2],[Hstar1_K1 Hstar1_K2],'o','Color',blue,'MarkerFaceColor','w', ...
    'MarkerSize',marker_size+1,'LineWidth',1.2);
plot([100*dK1 100*dK2],[Hstar2_K1 Hstar2_K2],'o','Color',orange,'MarkerFaceColor','w', ...
    'MarkerSize',marker_size+1,'LineWidth',1.2);
hold off;
xlabel('Decline in K (% per decade)','FontSize',label_size);
ylabel('Harvest threshold, H* (bears yr^{-1})','FontSize',label_size);
xlim([0 100*dK2]); ylim([0 125]);
legend('Lower-growth','Higher-growth', ...
    sprintf('Reference harvest: %.1f',H_recent), ...
    sprintf('Management reference: %.0f',H_management), ...
    'Location','southwest','FontSize',legend_size,'Box','off');
text(100*dK1,121,sprintf('%.0f%%',100*dK1),'HorizontalAlignment','center','VerticalAlignment','top', ...
    'FontSize',8,'Color',midgray);
text(100*dK2,121,sprintf('%.0f%%',100*dK2),'HorizontalAlignment','right','VerticalAlignment','top', ...
    'FontSize',8,'Color',midgray);

set(gca,'FontName',font_name,'FontSize',font_size,'LineWidth',0.8, ...
    'TickDir','out','Box','off','XColor',charcoal,'YColor',charcoal);
grid on; ax = gca; ax.GridColor = lightgray; ax.GridAlpha = 0.42;

drawnow;
savefig(f,fullfile(outdir,'Figure5_harvest_threshold.fig'));
set(f,'PaperPositionMode','auto');
print(f,fullfile(outdir,'Figure5_harvest_threshold.png'),'-dpng','-r600');



%% Supplementary Figure S1

N_ipm = 2937;
N_ipm_low = 1552;
N_ipm_high = 5944;

N_aerial1 = 3435;
N_aerial1_low = 2300;
N_aerial1_high = 5131;

N_aerial2 = 4196;
N_aerial2_low = 2807;
N_aerial2_high = 6273;

N_aerial3 = 5444;
N_aerial3_low = 3636;
N_aerial3_high = 8152;

estimates = [N_ipm N_aerial1 N_aerial2 N_aerial3];
lower_error = [N_ipm-N_ipm_low N_aerial1-N_aerial1_low ...
    N_aerial2-N_aerial2_low N_aerial3-N_aerial3_low];
upper_error = [N_ipm_high-N_ipm N_aerial1_high-N_aerial1 ...
    N_aerial2_high-N_aerial2 N_aerial3_high-N_aerial3];

f = figure('Color','w','Units','centimeters','Position',[2 2 14.5 8.2]);
set(f,'Renderer','painters');

hold on;
errorbar(1,estimates(1),lower_error(1),upper_error(1),'o', ...
    'Color',blue,'MarkerFaceColor',blue,'MarkerSize',7,'LineWidth',1.4,'CapSize',8);
for q = 2:4
    errorbar(q,estimates(q),lower_error(q),upper_error(q),'o', ...
        'Color',midgray,'MarkerFaceColor','w','MarkerSize',7,'LineWidth',1.3,'CapSize',8);
end
hold off;
ylabel('Population estimate','FontSize',label_size);
set(gca,'XTick',1:4,'XTickLabel',{'IPM','Aerial g(0)=1.0','Aerial g(0)=0.8','Aerial g(0)=0.6'});
xlim([0.5 4.5]); ylim([1000 8500]);
set(gca,'FontName',font_name,'FontSize',font_size,'LineWidth',0.8, ...
    'TickDir','out','Box','off','XColor',charcoal,'YColor',charcoal);
grid on; ax = gca; ax.GridColor = lightgray; ax.GridAlpha = 0.42;

drawnow;
savefig(f,fullfile(outdir,'Supplementary_Figure_S1_population_size_comparison.fig'));
set(f,'PaperPositionMode','auto');
print(f,fullfile(outdir,'Supplementary_Figure_S1_population_size_comparison.png'),'-dpng','-r600');


%% Figure 6

N_bounds = [N0_low N0_high];
r_bounds1 = [rMNPL1_low rMNPL1_high];
r_bounds2 = [rMNPL2_low rMNPL2_high];
x_bounds1 = [xMNPL1_low xMNPL1_high];
x_bounds2 = [xMNPL2_low xMNPL2_high];
dK_bounds = [dK1 dK2];
H_bounds = [H_recent H_high];

envelope1 = [];
envelope2 = [];

% Lower-growth scenario
for a = 1:2

    Nstart = N_bounds(a);

    for b = 1:2

        rMNPL = r_bounds1(b);

        for c = 1:2

            xMNPL = x_bounds1(c);
            fun = @(theta) (1/(theta+1))^(1/theta)-xMNPL;
            theta = fzero(fun,[0.01 50]);

            if rMNPL == 0
                r = 0;
            else
                r = rMNPL/(1-xMNPL^theta);
            end

            K0 = Nstart/xMNPL;

            for d = 1:2

                dK = dK_bounds(d);

                for e = 1:2

                    H = H_bounds(e);
                    N = zeros(size(tspan));
                    N(1) = Nstart;

                    for i = 2:length(tspan)

                        t = tspan(i-1);
                        odefun = @(tt,y) r*max(y,0)*(1-(max(y,0)/(K0*(1-dK)^(tt/10)))^theta)-H;

                        k1 = odefun(t,N(i-1));
                        k2 = odefun(t+dt/2,N(i-1)+dt*k1/2);
                        k3 = odefun(t+dt/2,N(i-1)+dt*k2/2);
                        k4 = odefun(t+dt,N(i-1)+dt*k3);

                        y = N(i-1)+dt*(k1+2*k2+2*k3+k4)/6;
                        if y < 0
                            y = 0;
                        end
                        N(i) = y;

                    end

                    envelope1 = [envelope1; N];

                end

            end

        end

    end

end

% Higher-growth scenario
for a = 1:2

    Nstart = N_bounds(a);

    for b = 1:2

        rMNPL = r_bounds2(b);

        for c = 1:2

            xMNPL = x_bounds2(c);
            fun = @(theta) (1/(theta+1))^(1/theta)-xMNPL;
            theta = fzero(fun,[0.01 50]);

            if rMNPL == 0
                r = 0;
            else
                r = rMNPL/(1-xMNPL^theta);
            end

            K0 = Nstart/xMNPL;

            for d = 1:2

                dK = dK_bounds(d);

                for e = 1:2

                    H = H_bounds(e);
                    N = zeros(size(tspan));
                    N(1) = Nstart;

                    for i = 2:length(tspan)

                        t = tspan(i-1);
                        odefun = @(tt,y) r*max(y,0)*(1-(max(y,0)/(K0*(1-dK)^(tt/10)))^theta)-H;

                        k1 = odefun(t,N(i-1));
                        k2 = odefun(t+dt/2,N(i-1)+dt*k1/2);
                        k3 = odefun(t+dt/2,N(i-1)+dt*k2/2);
                        k4 = odefun(t+dt,N(i-1)+dt*k3);

                        y = N(i-1)+dt*(k1+2*k2+2*k3+k4)/6;
                        if y < 0
                            y = 0;
                        end
                        N(i) = y;

                    end

                    envelope2 = [envelope2; N];

                end

            end

        end

    end

end

low_envelope1 = min(envelope1,[],1);
high_envelope1 = max(envelope1,[],1);
low_envelope2 = min(envelope2,[],1);
high_envelope2 = max(envelope2,[],1);

central1 = N_habitat1(2,:);
central2 = N_habitat2(2,:);

f = figure('Color','w','Units','centimeters','Position',[2 2 18 8.4]);
set(f,'Renderer','painters');

subplot(1,2,1);
fill([tspan fliplr(tspan)],[low_envelope1 fliplr(high_envelope1)], ...
    verylight,'EdgeColor','none'); hold on;
plot(tspan,central1,'Color',blue,'LineWidth',line_width2);
plot([tspan(1) tspan(end)],[N0 N0],'--','Color',midgray,'LineWidth',1.0);
hold off;
title('Lower-growth','FontSize',title_size,'FontWeight','bold');
xlabel('Years','FontSize',label_size);
ylabel('Population size','FontSize',label_size);
xlim([tspan(1) tspan(end)]); ylim([0 9000]);
text(0.02,0.96,'(a)','Units','normalized','FontWeight','bold','FontSize',title_size, ...
    'VerticalAlignment','top');

subplot(1,2,2);
fill([tspan fliplr(tspan)],[low_envelope2 fliplr(high_envelope2)], ...
    verylight,'EdgeColor','none'); hold on;
plot(tspan,central2,'Color',orange,'LineWidth',line_width2);
plot([tspan(1) tspan(end)],[N0 N0],'--','Color',midgray,'LineWidth',1.0);
hold off;
title('Higher-growth','FontSize',title_size,'FontWeight','bold');
xlabel('Years','FontSize',label_size);
ylabel('Population size','FontSize',label_size);
xlim([tspan(1) tspan(end)]); ylim([0 9000]);
text(0.02,0.96,'(b)','Units','normalized','FontWeight','bold','FontSize',title_size, ...
    'VerticalAlignment','top');

axs = findall(f,'Type','axes');
set(axs,'FontName',font_name,'FontSize',font_size,'LineWidth',0.8, ...
    'TickDir','out','Box','off','XColor',charcoal,'YColor',charcoal);
for k = 1:length(axs)
    grid(axs(k),'on');
    axs(k).GridColor = lightgray;
    axs(k).GridAlpha = 0.38;
end

drawnow;
savefig(f,fullfile(outdir,'Figure6_stress_test_envelope.fig'));
set(f,'PaperPositionMode','auto');
print(f,fullfile(outdir,'Figure6_stress_test_envelope.png'),'-dpng','-r600');


%% Supplementary Figure S2

% Reference projections are H = 44.4 and dK = 3% per decade.
baseline1 = N_habitat1(2,end);
baseline2 = N_habitat2(2,end);

sensitivity1 = zeros(3,2);
sensitivity2 = zeros(3,2);

% Initial population size, lower-growth scenario
for q = 1:2

    Nstart = N_bounds(q);
    K0 = Nstart/xMNPL1;
    N = zeros(size(tspan));
    N(1) = Nstart;

    for i = 2:length(tspan)

        t = tspan(i-1);
        odefun = @(tt,y) r1*max(y,0)*(1-(max(y,0)/(K0*(1-dK1)^(tt/10)))^theta1)-H_recent;

        k1 = odefun(t,N(i-1));
        k2 = odefun(t+dt/2,N(i-1)+dt*k1/2);
        k3 = odefun(t+dt/2,N(i-1)+dt*k2/2);
        k4 = odefun(t+dt,N(i-1)+dt*k3);

        y = N(i-1)+dt*(k1+2*k2+2*k3+k4)/6;
        if y < 0
            y = 0;
        end
        N(i) = y;

    end

    sensitivity1(1,q) = 100*(N(end)-baseline1)/baseline1;

end

% Growth at MNPL, lower-growth scenario
for q = 1:2

    rMNPL = r_bounds1(q);
    r = rMNPL/(1-xMNPL1^theta1);
    N = zeros(size(tspan));
    N(1) = N0;

    for i = 2:length(tspan)

        t = tspan(i-1);
        odefun = @(tt,y) r*max(y,0)*(1-(max(y,0)/(K1*(1-dK1)^(tt/10)))^theta1)-H_recent;

        k1 = odefun(t,N(i-1));
        k2 = odefun(t+dt/2,N(i-1)+dt*k1/2);
        k3 = odefun(t+dt/2,N(i-1)+dt*k2/2);
        k4 = odefun(t+dt,N(i-1)+dt*k3);

        y = N(i-1)+dt*(k1+2*k2+2*k3+k4)/6;
        if y < 0
            y = 0;
        end
        N(i) = y;

    end

    sensitivity1(2,q) = 100*(N(end)-baseline1)/baseline1;

end

% MNPL relative density, lower-growth scenario
for q = 1:2

    xMNPL = x_bounds1(q);
    fun = @(theta) (1/(theta+1))^(1/theta)-xMNPL;
    theta = fzero(fun,[0.01 50]);
    r = rMNPL1/(1-xMNPL^theta);
    K0 = N0/xMNPL;
    N = zeros(size(tspan));
    N(1) = N0;

    for i = 2:length(tspan)

        t = tspan(i-1);
        odefun = @(tt,y) r*max(y,0)*(1-(max(y,0)/(K0*(1-dK1)^(tt/10)))^theta)-H_recent;

        k1 = odefun(t,N(i-1));
        k2 = odefun(t+dt/2,N(i-1)+dt*k1/2);
        k3 = odefun(t+dt/2,N(i-1)+dt*k2/2);
        k4 = odefun(t+dt,N(i-1)+dt*k3);

        y = N(i-1)+dt*(k1+2*k2+2*k3+k4)/6;
        if y < 0
            y = 0;
        end
        N(i) = y;

    end

    sensitivity1(3,q) = 100*(N(end)-baseline1)/baseline1;

end

% Initial population size, higher-growth scenario
for q = 1:2

    Nstart = N_bounds(q);
    K0 = Nstart/xMNPL2;
    N = zeros(size(tspan));
    N(1) = Nstart;

    for i = 2:length(tspan)

        t = tspan(i-1);
        odefun = @(tt,y) r2*max(y,0)*(1-(max(y,0)/(K0*(1-dK1)^(tt/10)))^theta2)-H_recent;

        k1 = odefun(t,N(i-1));
        k2 = odefun(t+dt/2,N(i-1)+dt*k1/2);
        k3 = odefun(t+dt/2,N(i-1)+dt*k2/2);
        k4 = odefun(t+dt,N(i-1)+dt*k3);

        y = N(i-1)+dt*(k1+2*k2+2*k3+k4)/6;
        if y < 0
            y = 0;
        end
        N(i) = y;

    end

    sensitivity2(1,q) = 100*(N(end)-baseline2)/baseline2;

end

% Growth at MNPL, higher-growth scenario
for q = 1:2

    rMNPL = r_bounds2(q);
    r = rMNPL/(1-xMNPL2^theta2);
    N = zeros(size(tspan));
    N(1) = N0;

    for i = 2:length(tspan)

        t = tspan(i-1);
        odefun = @(tt,y) r*max(y,0)*(1-(max(y,0)/(K2*(1-dK1)^(tt/10)))^theta2)-H_recent;

        k1 = odefun(t,N(i-1));
        k2 = odefun(t+dt/2,N(i-1)+dt*k1/2);
        k3 = odefun(t+dt/2,N(i-1)+dt*k2/2);
        k4 = odefun(t+dt,N(i-1)+dt*k3);

        y = N(i-1)+dt*(k1+2*k2+2*k3+k4)/6;
        if y < 0
            y = 0;
        end
        N(i) = y;

    end

    sensitivity2(2,q) = 100*(N(end)-baseline2)/baseline2;

end

% MNPL relative density, higher-growth scenario
for q = 1:2

    xMNPL = x_bounds2(q);
    fun = @(theta) (1/(theta+1))^(1/theta)-xMNPL;
    theta = fzero(fun,[0.01 50]);
    r = rMNPL2/(1-xMNPL^theta);
    K0 = N0/xMNPL;
    N = zeros(size(tspan));
    N(1) = N0;

    for i = 2:length(tspan)

        t = tspan(i-1);
        odefun = @(tt,y) r*max(y,0)*(1-(max(y,0)/(K0*(1-dK1)^(tt/10)))^theta)-H_recent;

        k1 = odefun(t,N(i-1));
        k2 = odefun(t+dt/2,N(i-1)+dt*k1/2);
        k3 = odefun(t+dt/2,N(i-1)+dt*k2/2);
        k4 = odefun(t+dt,N(i-1)+dt*k3);

        y = N(i-1)+dt*(k1+2*k2+2*k3+k4)/6;
        if y < 0
            y = 0;
        end
        N(i) = y;

    end

    sensitivity2(3,q) = 100*(N(end)-baseline2)/baseline2;

end

f = figure('Color','w','Units','centimeters','Position',[2 2 18 8.4]);
set(f,'Renderer','painters');

subplot(1,2,1); hold on;
for q = 1:3
    plot(sensitivity1(q,:),[q q],'Color',blue,'LineWidth',line_width2);
    plot(sensitivity1(q,1),q,'o','Color',blue,'MarkerFaceColor','w','MarkerSize',marker_size+1,'LineWidth',1.2);
    plot(sensitivity1(q,2),q,'s','Color',blue,'MarkerFaceColor',blue,'MarkerSize',marker_size,'LineWidth',1.1);
end
plot([0 0],[0.5 3.5],'--','Color',midgray,'LineWidth',1.0);
hold off;
title('Lower-growth','FontSize',title_size,'FontWeight','bold');
xlabel('Change in N(35) from reference (%)','FontSize',label_size);
set(gca,'YTick',1:3,'YTickLabel',{'Initial population size','Growth at MNPL','MNPL relative density'}, ...
    'YDir','reverse');
ylim([0.5 3.5]);
text(0.02,0.96,'(a)','Units','normalized','FontWeight','bold','FontSize',title_size, ...
    'VerticalAlignment','top');

subplot(1,2,2); hold on;
for q = 1:3
    plot(sensitivity2(q,:),[q q],'Color',orange,'LineWidth',line_width2);
    plot(sensitivity2(q,1),q,'o','Color',orange,'MarkerFaceColor','w','MarkerSize',marker_size+1,'LineWidth',1.2);
    plot(sensitivity2(q,2),q,'s','Color',orange,'MarkerFaceColor',orange,'MarkerSize',marker_size,'LineWidth',1.1);
end
plot([0 0],[0.5 3.5],'--','Color',midgray,'LineWidth',1.0);
hold off;
title('Higher-growth','FontSize',title_size,'FontWeight','bold');
xlabel('Change in N(35) from reference (%)','FontSize',label_size);
set(gca,'YTick',1:3,'YTickLabel',{'Initial population size','Growth at MNPL','MNPL relative density'}, ...
    'YDir','reverse');
ylim([0.5 3.5]);
text(0.02,0.96,'(b)','Units','normalized','FontWeight','bold','FontSize',title_size, ...
    'VerticalAlignment','top');

axs = findall(f,'Type','axes');
set(axs,'FontName',font_name,'FontSize',font_size,'LineWidth',0.8, ...
    'TickDir','out','Box','off','XColor',charcoal,'YColor',charcoal);
for k = 1:length(axs)
    grid(axs(k),'on');
    axs(k).GridColor = lightgray;
    axs(k).GridAlpha = 0.42;
end

drawnow;
savefig(f,fullfile(outdir,'Supplementary_Figure_S2_OAT_sensitivity.fig'));
set(f,'PaperPositionMode','auto');
print(f,fullfile(outdir,'Supplementary_Figure_S2_OAT_sensitivity.png'),'-dpng','-r600');


%% Supplementary Figure S3

initial_density = 0.50:0.01:0.94;

N_density1_K1 = zeros(size(initial_density));
N_density1_K2 = zeros(size(initial_density));
N_density2_K1 = zeros(size(initial_density));
N_density2_K2 = zeros(size(initial_density));

H_density1_K1 = nan(size(initial_density));
H_density1_K2 = nan(size(initial_density));
H_density2_K1 = nan(size(initial_density));
H_density2_K2 = nan(size(initial_density));

% 35 year population for different starting relative densities
for q = 1:length(initial_density)

    x0 = initial_density(q);
    K0_1 = N0/x0;
    K0_2 = N0/x0;

    % Lower-growth, 3% decline
    dK = dK1;
    H = H_recent;
    N = zeros(size(tspan));
    N(1) = N0;

    for i = 2:length(tspan)
        t = tspan(i-1);
        odefun = @(tt,y) r1*max(y,0)*(1-(max(y,0)/(K0_1*(1-dK)^(tt/10)))^theta1)-H;
        k1 = odefun(t,N(i-1));
        k2 = odefun(t+dt/2,N(i-1)+dt*k1/2);
        k3 = odefun(t+dt/2,N(i-1)+dt*k2/2);
        k4 = odefun(t+dt,N(i-1)+dt*k3);
        y = N(i-1)+dt*(k1+2*k2+2*k3+k4)/6;
        if y < 0
            y = 0;
        end
        N(i) = y;
    end

    N_density1_K1(q) = N(end);

    % Lower-growth, 9% decline
    dK = dK2;
    H = H_recent;
    N = zeros(size(tspan));
    N(1) = N0;

    for i = 2:length(tspan)
        t = tspan(i-1);
        odefun = @(tt,y) r1*max(y,0)*(1-(max(y,0)/(K0_1*(1-dK)^(tt/10)))^theta1)-H;
        k1 = odefun(t,N(i-1));
        k2 = odefun(t+dt/2,N(i-1)+dt*k1/2);
        k3 = odefun(t+dt/2,N(i-1)+dt*k2/2);
        k4 = odefun(t+dt,N(i-1)+dt*k3);
        y = N(i-1)+dt*(k1+2*k2+2*k3+k4)/6;
        if y < 0
            y = 0;
        end
        N(i) = y;
    end

    N_density1_K2(q) = N(end);

    % Higher-growth, 3% decline
    dK = dK1;
    H = H_recent;
    N = zeros(size(tspan));
    N(1) = N0;

    for i = 2:length(tspan)
        t = tspan(i-1);
        odefun = @(tt,y) r2*max(y,0)*(1-(max(y,0)/(K0_2*(1-dK)^(tt/10)))^theta2)-H;
        k1 = odefun(t,N(i-1));
        k2 = odefun(t+dt/2,N(i-1)+dt*k1/2);
        k3 = odefun(t+dt/2,N(i-1)+dt*k2/2);
        k4 = odefun(t+dt,N(i-1)+dt*k3);
        y = N(i-1)+dt*(k1+2*k2+2*k3+k4)/6;
        if y < 0
            y = 0;
        end
        N(i) = y;
    end

    N_density2_K1(q) = N(end);

    % Higher-growth, 9% decline
    dK = dK2;
    H = H_recent;
    N = zeros(size(tspan));
    N(1) = N0;

    for i = 2:length(tspan)
        t = tspan(i-1);
        odefun = @(tt,y) r2*max(y,0)*(1-(max(y,0)/(K0_2*(1-dK)^(tt/10)))^theta2)-H;
        k1 = odefun(t,N(i-1));
        k2 = odefun(t+dt/2,N(i-1)+dt*k1/2);
        k3 = odefun(t+dt/2,N(i-1)+dt*k2/2);
        k4 = odefun(t+dt,N(i-1)+dt*k3);
        y = N(i-1)+dt*(k1+2*k2+2*k3+k4)/6;
        if y < 0
            y = 0;
        end
        N(i) = y;
    end

    N_density2_K2(q) = N(end);

end

% Thresholds for different starting relative densities.
H_density_max = 200;

for q = 1:length(initial_density)

    x0 = initial_density(q);
    K0_1 = N0/x0;
    K0_2 = N0/x0;

    for scenario = 1:4

        if scenario == 1
            r = r1;
            theta = theta1;
            K0 = K0_1;
            dK = dK1;
        elseif scenario == 2
            r = r1;
            theta = theta1;
            K0 = K0_1;
            dK = dK2;
        elseif scenario == 3
            r = r2;
            theta = theta2;
            K0 = K0_2;
            dK = dK1;
        else
            r = r2;
            theta = theta2;
            K0 = K0_2;
            dK = dK2;
        end

        Hlo = 0;
        Hhi = H_density_max;

        H = Hlo;
        N = zeros(size(tspan));
        N(1) = N0;

        for i = 2:length(tspan)
            t = tspan(i-1);
            odefun = @(tt,y) r*max(y,0)*(1-(max(y,0)/(K0*(1-dK)^(tt/10)))^theta)-H;
            k1 = odefun(t,N(i-1));
            k2 = odefun(t+dt/2,N(i-1)+dt*k1/2);
            k3 = odefun(t+dt/2,N(i-1)+dt*k2/2);
            k4 = odefun(t+dt,N(i-1)+dt*k3);
            y = N(i-1)+dt*(k1+2*k2+2*k3+k4)/6;
            if y < 0
                y = 0;
            end
            N(i) = y;
        end

        Flo = N(end)-N0;

        H = Hhi;
        N = zeros(size(tspan));
        N(1) = N0;

        for i = 2:length(tspan)
            t = tspan(i-1);
            odefun = @(tt,y) r*max(y,0)*(1-(max(y,0)/(K0*(1-dK)^(tt/10)))^theta)-H;
            k1 = odefun(t,N(i-1));
            k2 = odefun(t+dt/2,N(i-1)+dt*k1/2);
            k3 = odefun(t+dt/2,N(i-1)+dt*k2/2);
            k4 = odefun(t+dt,N(i-1)+dt*k3);
            y = N(i-1)+dt*(k1+2*k2+2*k3+k4)/6;
            if y < 0
                y = 0;
            end
            N(i) = y;
        end

        Fhi = N(end)-N0;
        Hresult = nan;

        if Flo >= 0 && Fhi <= 0

            while Hhi-Hlo > threshold_tolerance

                Hmid = (Hlo+Hhi)/2;
                H = Hmid;
                N = zeros(size(tspan));
                N(1) = N0;

                for i = 2:length(tspan)
                    t = tspan(i-1);
                    odefun = @(tt,y) r*max(y,0)*(1-(max(y,0)/(K0*(1-dK)^(tt/10)))^theta)-H;
                    k1 = odefun(t,N(i-1));
                    k2 = odefun(t+dt/2,N(i-1)+dt*k1/2);
                    k3 = odefun(t+dt/2,N(i-1)+dt*k2/2);
                    k4 = odefun(t+dt,N(i-1)+dt*k3);
                    y = N(i-1)+dt*(k1+2*k2+2*k3+k4)/6;
                    if y < 0
                        y = 0;
                    end
                    N(i) = y;
                end

                Fmid = N(end)-N0;

                if Fmid >= 0
                    Hlo = Hmid;
                else
                    Hhi = Hmid;
                end

            end

            Hresult = Hlo;

        end

        if scenario == 1
            H_density1_K1(q) = Hresult;
        elseif scenario == 2
            H_density1_K2(q) = Hresult;
        elseif scenario == 3
            H_density2_K1(q) = Hresult;
        else
            H_density2_K2(q) = Hresult;
        end

    end

end

f = figure('Color','w','Units','centimeters','Position',[2 2 18 8.4]);
set(f,'Renderer','painters');

subplot(1,2,1); hold on;
plot(initial_density,N_density1_K1,'Color',blue,'LineStyle','-','LineWidth',line_width2);
plot(initial_density,N_density1_K2,'Color',blue,'LineStyle','--','LineWidth',line_width2);
plot(initial_density,N_density2_K1,'Color',orange,'LineStyle','-','LineWidth',line_width2);
plot(initial_density,N_density2_K2,'Color',orange,'LineStyle','--','LineWidth',line_width2);
plot([0.50 0.94],[N0 N0],'--','Color',midgray,'LineWidth',1.0);
hold off;
xlabel('Initial relative density, N_0/K_0','FontSize',label_size);
ylabel('Population size after 35 years','FontSize',label_size);
legend(sprintf('Lower-growth, %.0f%%',100*dK1),sprintf('Lower-growth, %.0f%%',100*dK2), ...
    sprintf('Higher-growth, %.0f%%',100*dK1),sprintf('Higher-growth, %.0f%%',100*dK2), ...
    'Location','northeast','FontSize',legend_size,'Box','off');
xlim([initial_density(1) initial_density(end)]);
text(0.02,0.96,'(a)','Units','normalized','FontWeight','bold','FontSize',title_size, ...
    'VerticalAlignment','top');

subplot(1,2,2); hold on;
plot(initial_density,H_density1_K1,'Color',blue,'LineStyle','-','LineWidth',line_width2);
plot(initial_density,H_density1_K2,'Color',blue,'LineStyle','--','LineWidth',line_width2);
plot(initial_density,H_density2_K1,'Color',orange,'LineStyle','-','LineWidth',line_width2);
plot(initial_density,H_density2_K2,'Color',orange,'LineStyle','--','LineWidth',line_width2);
hold off;
xlabel('Initial relative density, N_0/K_0','FontSize',label_size);
ylabel('Harvest threshold, H* (bears yr^{-1})','FontSize',label_size);
legend(sprintf('Lower-growth, %.0f%%',100*dK1),sprintf('Lower-growth, %.0f%%',100*dK2), ...
    sprintf('Higher-growth, %.0f%%',100*dK1),sprintf('Higher-growth, %.0f%%',100*dK2), ...
    'Location','northeast','FontSize',legend_size,'Box','off');
xlim([initial_density(1) initial_density(end)]); ylim([0 150]);
text(0.02,0.96,'(b)','Units','normalized','FontWeight','bold','FontSize',title_size, ...
    'VerticalAlignment','top');

axs = findall(f,'Type','axes');
set(axs,'FontName',font_name,'FontSize',font_size,'LineWidth',0.8, ...
    'TickDir','out','Box','off','XColor',charcoal,'YColor',charcoal);
for k = 1:length(axs)
    grid(axs(k),'on');
    axs(k).GridColor = lightgray;
    axs(k).GridAlpha = 0.42;
end

drawnow;
savefig(f,fullfile(outdir,'Supplementary_Figure_S3_initial_density.fig'));
set(f,'PaperPositionMode','auto');
print(f,fullfile(outdir,'Supplementary_Figure_S3_initial_density.png'),'-dpng','-r600');


%% Numerical convergence

dt_check = 0.025;
tspan_check = 0:dt_check:35;

% Lower-growth scenario, contemporary harvesting and 3% K decline
N_check1 = zeros(size(tspan_check));
N_check1(1) = N0;

for i = 2:length(tspan_check)

    t = tspan_check(i-1);
    odefun = @(tt,y) r1*max(y,0)*(1-(max(y,0)/(K1*(1-dK1)^(tt/10)))^theta1)-H_recent;

    k1 = odefun(t,N_check1(i-1));
    k2 = odefun(t+dt_check/2,N_check1(i-1)+dt_check*k1/2);
    k3 = odefun(t+dt_check/2,N_check1(i-1)+dt_check*k2/2);
    k4 = odefun(t+dt_check,N_check1(i-1)+dt_check*k3);

    y = N_check1(i-1)+dt_check*(k1+2*k2+2*k3+k4)/6;
    if y < 0
        y = 0;
    end
    N_check1(i) = y;

end

% Higher-growth scenario, contemporary harvesting and 3% K decline
N_check2 = zeros(size(tspan_check));
N_check2(1) = N0;

for i = 2:length(tspan_check)

    t = tspan_check(i-1);
    odefun = @(tt,y) r2*max(y,0)*(1-(max(y,0)/(K2*(1-dK1)^(tt/10)))^theta2)-H_recent;

    k1 = odefun(t,N_check2(i-1));
    k2 = odefun(t+dt_check/2,N_check2(i-1)+dt_check*k1/2);
    k3 = odefun(t+dt_check/2,N_check2(i-1)+dt_check*k2/2);
    k4 = odefun(t+dt_check,N_check2(i-1)+dt_check*k3);

    y = N_check2(i-1)+dt_check*(k1+2*k2+2*k3+k4)/6;
    if y < 0
        y = 0;
    end
    N_check2(i) = y;

end

difference1 = abs(N_habitat1(2,end)-N_check1(end));
difference2 = abs(N_habitat2(2,end)-N_check2(end));
relative_difference1 = 100*difference1/N_check1(end);
relative_difference2 = 100*difference2/N_check2(end);


%% Summary of results

fprintf('\n--------------------------------------------------\n');
fprintf('35-year population projections\n');
fprintf('--------------------------------------------------\n');
fprintf('Lower-growth, 3%% K decline: %.3f bears (%.1f%% of N0)\n', ...
    N_habitat1(2,end),100*N_habitat1(2,end)/N0);
fprintf('Lower-growth, 9%% K decline: %.3f bears (%.1f%% of N0)\n', ...
    N_habitat1(3,end),100*N_habitat1(3,end)/N0);
fprintf('Higher-growth, 3%% K decline: %.3f bears (%.1f%% of N0)\n', ...
    N_habitat2(2,end),100*N_habitat2(2,end)/N0);
fprintf('Higher-growth, 9%% K decline: %.3f bears (%.1f%% of N0)\n', ...
    N_habitat2(3,end),100*N_habitat2(3,end)/N0);

fprintf('\nInitial net population growth under historical harvesting of 56 bears yr^-1\n');
fprintf('Published estimate: %.3f yr^-1 (95%% CRI %.3f to %.3f)\n', ...
    g_observed,g_observed_low,g_observed_high);
fprintf('Lower-growth model: %.3f yr^-1\n',g_model1);
fprintf('Higher-growth model: %.3f yr^-1\n',g_model2);


fprintf('\nAnalytical maximum-production ceiling\n');
fprintf('Initial Pmax: %.1f bears yr^-1 (lower) and %.1f bears yr^-1 (higher)\n', ...
    Pmax1_0,Pmax2_0);
fprintf('Pmax at 35 years, lower-growth: %.1f (3%%) and %.1f (9%%) bears yr^-1\n', ...
    Pmax1_K1_35,Pmax1_K2_35);
fprintf('Pmax at 35 years, higher-growth: %.1f (3%%) and %.1f (9%%) bears yr^-1\n', ...
    Pmax2_K1_35,Pmax2_K2_35);
fprintf('Reference harvest crosses lower-growth 9%% Pmax at %.1f years\n',tcross1_K2);
fprintf('Other reference crossings occur after the 35-year projection horizon\n');

fprintf('\nHarvest threshold for N(35) = N0\n');
fprintf('Lower-growth, 3%% K decline: %.1f bears yr^-1\n',Hstar1_K1);
fprintf('Higher-growth, 3%% K decline: %.1f bears yr^-1\n',Hstar2_K1);
fprintf('Lower-growth, 9%% K decline: %.1f bears yr^-1\n',Hstar1_K2);
fprintf('Higher-growth, 9%% K decline: %.1f bears yr^-1\n',Hstar2_K2);

fprintf('\nStructural sensitivity to initial relative density\n');
fprintf('N0/K0 = 0.50 to 0.94\n');
fprintf('Lower-growth, 3%%: %.0f to %.0f bears\n', ...
    N_density1_K1(1),N_density1_K1(end));
fprintf('Lower-growth, 9%%: %.0f to %.0f bears\n', ...
    N_density1_K2(1),N_density1_K2(end));
fprintf('Higher-growth, 3%%: %.0f to %.0f bears\n', ...
    N_density2_K1(1),N_density2_K1(end));
fprintf('Higher-growth, 9%%: %.0f to %.0f bears\n', ...
    N_density2_K2(1),N_density2_K2(end));

fprintf('\nNumerical convergence\n');
fprintf('Lower-growth: dt=0.05 %.6f, dt=0.025 %.6f, relative difference %.12f %%\n', ...
    N_habitat1(2,end),N_check1(end),relative_difference1);
fprintf('Higher-growth: dt=0.05 %.6f, dt=0.025 %.6f, relative difference %.12f %%\n', ...
    N_habitat2(2,end),N_check2(end),relative_difference2);
