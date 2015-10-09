height=0.74;
width=0.4;

fig=figure('Color',[1 1 1]);
clf;

set(fig, 'Position', [100 100 800 300])

%General parameters
nmtds=3;


names1=cellstr(char('TRSK-HCT-HR95', ...
    'MODF-HCM-HR95'));

format short e
marcvec=cellstr(char('-*r','-og','-sb','-^c','-+m','-dy','-vk','-pk'));

marksize=1;
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

tfinal=1200;

%----------------
%1st graph
%-------------
subplot1=subplot(1,2,1);
set(subplot1, 'Position',[0.08 0.2 width height]);
x=ds.timedys(1:tfinal);
y=ds.trskerrormax_h(1:tfinal);
i=1
plot(x,y ,[linestyles{i} Markers(i)],'LineWidth',linewidth, ...
    'MarkerSize', marksize, 'DisplayName', char(names1(i)), 'Color', colors(i,:))
hold on
x=ds.timedys(1:tfinal);
y=ds.newerrormax_h(1:tfinal);
i=2
plot(x,y ,[linestyles{i} Markers(i)],'LineWidth',linewidth, ...
    'MarkerSize', marksize, 'DisplayName', char(names1(i)), 'Color', colors(i,:))

ylabel('Max Error','fontsize',12)
xlabel('Days','fontsize',12)
%title('H','fontsize',12 )
box off;

%----------------
%2st graph
%-------------
subplot2=subplot(1,2,2);
set(subplot2, 'Position',[0.58 0.2 width height]);

x=ds.timedys(1:tfinal);
y=ds.trskerror2_h(1:tfinal);
i=1
plot(x,y ,[linestyles{i} Markers(i)],'LineWidth',linewidth, ...
    'MarkerSize', marksize, 'DisplayName', char(names1(i)), 'Color', colors(i,:))
hold on
x=ds.timedys(1:tfinal);
y=ds.newerror2_h(1:tfinal);
i=2
plot(x,y ,[linestyles{i} Markers(i)],'LineWidth',linewidth, ...
    'MarkerSize', marksize, 'DisplayName', char(names1(i)), 'Color', colors(i,:))



ylabel('RMS Error','fontsize',12)
xlabel('Days','fontsize',12)
%title(' RMS Errors','fontsize',12 )
box off;

leg2=legend(subplot2, 'show'); 
set(leg2, 'Position',[0.15 0.02 0.8 0.03],'Orientation','horizontal');
%set(leg2, 'Position',[0.84 0.45 0.13 0.13]);
legend boxoff;
