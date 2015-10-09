height=0.74;
width=0.4;

fig=figure('Color',[1 1 1]);
clf;

set(fig, 'Position', [100 100 800 300])

%General parameters
nmtds=3;

grids=cellstr(char('icos_pol_nopt_', ...
    'icos_pol_scvt_h1_',...
    'HR95JT_00'));

names1=cellstr(char(...
    'TRSK-HCT-HR95', ...
    'MODF-HCM-HR95' ...
     ));
%'TRSK-HCT-SCVT',...
%'MODF-HCM-HR95*'

format short e
marcvec=cellstr(char('-*r','-og','-sb','-^c','-+m','-dy','-vk','-pk'));

marksize=0.0001;
linewidth=1.2;
linestyles = cellstr(char('-','--','-.','-','-','--','-.',':','-',':','-',':',...
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

tfinal=600;

%----------------
%1st graph
%-------------
subplot1=subplot(1,2,1);
set(subplot1, 'Position',[0.08 0.16 width height]);
% x=tc2scvt.timedys(2:tfinal);
% y=tc2scvt.Kenergy(2:tfinal);
% i=1
% plot(x,y ,[linestyles{i} Markers(i)],'LineWidth',linewidth, ...
%     'MarkerSize', marksize, 'DisplayName', char(names1(i)), 'Color', colors(i,:))
% hold on

x=tc2trskHR957.timedys(2:tfinal);
y=tc2trskHR957.Kenergy(2:tfinal);
i=1
plot(x,y ,[linestyles{i} Markers(i)],'LineWidth',linewidth, ...
    'MarkerSize', marksize, 'DisplayName', char(names1(i)), 'Color', colors(i,:))
hold on 
x=tc2new7.timedys(2:tfinal);
y=tc2new7.Kenergy(2:tfinal);
i=2
plot(x,y ,[linestyles{i} Markers(i)],'LineWidth',linewidth, ...
    'MarkerSize', marksize, 'DisplayName', char(names1(i)), 'Color', colors(i,:))

%Proper calculated kinetic energy of the new method using perots method
% x=tc2new.timedys(2:tfinal);
% y=tc2new_kenergy.Kenergy(2:tfinal);
% i=4
% plot(x,y ,[linestyles{i} Markers(i)],'LineWidth',linewidth, ...
%     'MarkerSize', marksize, 'DisplayName', char(names1(i)), 'Color', colors(i,:))
% 
ylabel('Max Error','fontsize',12)
xlabel('Days','fontsize',12)
title('Kinetic Energy Variation','fontsize',12 )
box off;

%----------------
%2st graph
%-------------
subplot2=subplot(1,2,2);
set(subplot2, 'Position',[0.54 0.16 width height]);

% x=tc2scvt.timedys(2:tfinal);
% y=tc2scvt.Tenergy(2:tfinal);
% %y=abs(y);
% i=1
% plot(x,y ,[linestyles{i} Markers(i)],'LineWidth',linewidth, ...
%     'MarkerSize', marksize, 'DisplayName', char(names1(i)), 'Color', colors(i,:))

x=tc2trskHR957.timedys(2:tfinal);
y=tc2trskHR957.Tenergy(2:tfinal);
y=abs(y);
i=1
semilogy(x,y ,[linestyles{i} Markers(i)],'LineWidth',linewidth, ...
    'MarkerSize', marksize, 'DisplayName', char(names1(i)), 'Color', colors(i,:))
hold on


x=tc2new7.timedys(2:tfinal);
y=tc2new7.Tenergy(2:tfinal);
%y=abs(y);
i=2
semilogy(x,y ,[linestyles{i} Markers(i)],'LineWidth',linewidth, ...
    'MarkerSize', marksize, 'DisplayName', char(names1(i)), 'Color', colors(i,:))

%x=tc2new.timedys(2:tfinal);
%y=tc2new_kenergy.Tenergy(2:tfinal);
%y=abs(y);
% i=4
% plot(x,y ,[linestyles{i} Markers(i)],'LineWidth',linewidth, ...
%     'MarkerSize', marksize, 'DisplayName', char(names1(i)), 'Color', colors(i,:))

%ylabel('Max Error','fontsize',12)
xlabel('Days','fontsize',12)
title(' Total Energy Variation','fontsize',12 )
box off;

leg2=legend(subplot2, 'show'); 
%set(leg2, 'Position',[0.84 0.45 0.13 0.13]);
set(leg2, 'Position',[0.1 0.00 0.8 0.03],'Orientation','horizontal');
legend boxoff;
