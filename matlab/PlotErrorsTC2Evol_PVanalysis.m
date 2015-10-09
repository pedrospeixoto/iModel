height=0.35;
width=0.6;

fig=figure('Color',[1 1 1]);
clf;

set(fig, 'Position', [100 100 800 500])

%General parameters
nmtds=3;


% 
% names1=cellstr(char('TRSK-SCVT-ORGPV',...
%     'TRSK-HR95-ORGPV', ...
%     'TRSK-SCVT-CLUST', ...
%     'TRSK-HR95-CLUST', ...
%     'MODF-HR95-ORGPV', ...
%     'MODF-HR95-CLUST'));

names1=cellstr(char(...
    'TRSK-HR95-ORGPV', ...
    'MODF-HR95-ORGPV', ...
    'TRSK-HR95-APVM', ...
    'MODF-HR95-APVM'));

format short e
marcvec=cellstr(char('-*r','-og','-sb','-^c','-+m','-dy','-vk','-pk'));

marksize=0.0001;
linewidth=1.2;
linestyles = cellstr(char('--','-',':','-.','-','--','-.',':','-',':','-',':',...
    '-.','--','-',':','-.','--','-',':','-.'));


Markers=['o','s','+','x','^','d','v','^','<','>','p','h','.',...
    '+','*','o','x','^','<','h','.','>','p','s','d','v',...
    'o','x','+','*','s','d','v','^','<','>','p','h','.'];

%colors(1:7, 1:3)=0;
colors=jet(7);
colors(1,:)=[204/256 0/256 0/256];
colors(2,:)=[0/256 153/256 0/256];
colors(3,:)=[0/256 0/256 153/256];
colors(4,:)=[162/256 0/256 162/256];
colors(5,:)=[0/256 162/256 162/256];
colors(6,:)=[162/256 162/256 0/256];
colors(7,:)=[192/256 192/256 192/256];

tfinal=11999;

%----------------
%1st graph
%-------------
subplot1=subplot(2,1,1);
set(subplot1, 'Position',[0.08 0.59 width height]);

%x=tc2trskSCVT.timedys(1:end);
%y=tc2trskSCVT.errormax_h(1:end);
%i=1
%semilogy(x,y ,[linestyles{i} Markers(i)],'LineWidth',linewidth, ...
%    'MarkerSize', marksize, 'DisplayName', char(names1(i)), 'Color', colors(i,:))
%hold on

x=tc2trskHR957.timedys(1:end);
y=tc2trskHR957.errormax_h(1:end);
i=1
semilogy(x,y ,[linestyles{i} Markers(i)],'LineWidth',2.0, ...
    'MarkerSize', marksize, 'DisplayName', char(names1(i)), 'Color', colors(i,:))
hold on

%x=tc2trskSCVTclust.timedys(1:end);
%y=tc2trskSCVTclust.errormax_h(1:end);
%i=3
%semilogy(x,y ,[linestyles{i} Markers(i)],'LineWidth',linewidth, ...
%    'MarkerSize', marksize, 'DisplayName', char(names1(i)), 'Color', colors(i,:))

% x=tc2trskHR95clust.timedys(1:end);
% y=tc2trskHR95clust.errormax_h(1:end);
% i=2
% semilogy(x,y ,[linestyles{i} Markers(i)],'LineWidth',linewidth, ...
%     'MarkerSize', marksize, 'DisplayName', char(names1(i)), 'Color', colors(i,:))

x=tc2new7.timedys(1:end);
y=tc2new7.errormax_h(1:end);
i=2
semilogy(x,y ,[linestyles{i} Markers(i)],'LineWidth',linewidth, ...
    'MarkerSize', marksize, 'DisplayName', char(names1(i)), 'Color', colors(i,:))

x=tc2trskHR957apvm.timedys(1:end);
y=tc2trskHR957apvm.errormax_h(1:end);
i=3
semilogy(x,y ,[linestyles{i} Markers(i)],'LineWidth',linewidth, ...
    'MarkerSize', marksize, 'DisplayName', char(names1(i)), 'Color', colors(i,:))

x=tc2new7apvm.timedys(1:end);
y=tc2new7apvm.errormax_h(1:end);
i=4
semilogy(x,y ,[linestyles{i} Markers(i)],'LineWidth',linewidth, ...
    'MarkerSize', marksize, 'DisplayName', char(names1(i)), 'Color', colors(i,:))

axis([0 200 0.000001 1])
set(gca, 'fontsize', 12)

ylabel('Max Error','fontsize',12)
%xlabel('Days','fontsize',12)
title('H','fontsize',12 )
box off;

%----------------
%2st graph
%-------------
subplot2=subplot(2,1,2);
set(subplot2, 'Position',[0.08 0.09 width height]);

%x=tc2trskSCVT.timedys(1:end);
%y=tc2trskSCVT.errormax_pv(1:end);
%i=1
%semilogy(x,y ,[linestyles{i} Markers(i)],'LineWidth',linewidth, ...
%    'MarkerSize', marksize, 'DisplayName', char(names1(i)), 'Color', colors(i,:))
%hold on

x=tc2trskHR957.timedys(1:end);
y=tc2trskHR957.errormax_pv(1:end);
i=1
semilogy(x,y ,[linestyles{i} Markers(i)],'LineWidth',linewidth, ...
    'MarkerSize', marksize, 'DisplayName', char(names1(i)), 'Color', colors(i,:))
hold on

x=tc2new7.timedys(1:end);
y=tc2new7.errormax_pv(1:end);
i=2
semilogy(x,y ,[linestyles{i} Markers(i)],'LineWidth',linewidth, ...
    'MarkerSize', marksize, 'DisplayName', char(names1(i)), 'Color', colors(i,:))

%x=tc2trskSCVTclust.timedys(1:end);
%y=tc2trskSCVTclust.errormax_pv(1:end);
%i=3
%semilogy(x,y ,[linestyles{i} Markers(i)],'LineWidth',linewidth, ...
%    'MarkerSize', marksize, 'DisplayName', char(names1(i)), 'Color', colors(i,:))


% x=tc2trskHR95clust.timedys(1:end);
% y=tc2trskHR95clust.errormax_pv(1:end);
% i=2
% semilogy(x,y ,[linestyles{i} Markers(i)],'LineWidth',linewidth, ...
%     'MarkerSize', marksize, 'DisplayName', char(names1(i)), 'Color', colors(i,:))

x=tc2trskHR957apvm.timedys(1:end);
y=tc2trskHR957apvm.errormax_pv(1:end);
i=3
semilogy(x,y ,[linestyles{i} Markers(i)],'LineWidth',linewidth, ...
    'MarkerSize', marksize, 'DisplayName', char(names1(i)), 'Color', colors(i,:))

x=tc2new7apvm.timedys(1:end);
y=tc2new7apvm.errormax_pv(1:end);
i=4
semilogy(x,y ,[linestyles{i} Markers(i)],'LineWidth',linewidth, ...
    'MarkerSize', marksize, 'DisplayName', char(names1(i)), 'Color', colors(i,:))


ylabel('Max Error','fontsize',12)
xlabel('Days','fontsize',12)
title(' PV','fontsize',12 )
box off;

axis([0 200 0.00001 100])
set(gca, 'fontsize', 12)

leg2=legend(subplot2, 'show'); 
set(leg2, 'Position',[0.78 0.45 0.13 0.13],'fontsize',11);
legend boxoff;
