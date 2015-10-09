%% read_data %call read_data.m
height=0.80;
width=0.3;

fig=figure('Color',[1 1 1]);
clf;

set(fig, 'Position', [100 100 800 300])

%General parameters
nmtds=3;

grids=cellstr(char('icos_pol_nopt_', ...
    'icos_pol_scvt_h1_',...
    'HR95JT_00'));

names1=cellstr(char('ICOS',...
    'SCVT', ...
    'HR95'));


format short e
marcvec=cellstr(char('-*r','-og','-sb','-^c','-+m','-dy','-vk','-pk'));

marksize=6;
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


%----------------
%1st graph
%-------------
subplot1=subplot(1,2,1);
set(subplot1, 'Position',[0.08 0.12 width height]);

variable='errorinf'

%Reference lines
ref1=0.0025;
ref2=0.04;

hold off;
i=1;
yds=ds(strcmp(ds.grid,grids(i)) , {'Glevel',variable});
x=yds.Glevel(1:end);
y=yds(1:end,2); 
y

semilogy(x,y ,[linestyles{i} Markers(i)],'LineWidth',linewidth, ...
    'MarkerSize', marksize, 'DisplayName', char(names1(i)), 'Color', colors(i,:))
hold on

for i=2:3
yds=ds(strcmp(ds.grid,grids(i)) , {'Glevel',variable});
    x=yds.Glevel(1:end);
    y=yds(1:end,2); 
    y
semilogy(x, y,[linestyles{i} Markers(i)],'LineWidth',linewidth, 'MarkerSize', marksize, 'DisplayName', char(names1(i)), 'Color', colors(i,:))
end
n=size(y);
%Reference lines
ref(1,1)=ref1;
ref(1,2)=ref2;
xref(1)=1;
for j=2:size(y)
    xref(j)=j;
    ref(j, 1)=ref(j-1, 1)/2;
    ref(j, 2)=ref(j-1, 2)/4;
end
semilogy(xref(2:4),ref(2:4,1),'-', 'Color', [152/256 152/256 152/256],'LineWidth',1.0, 'DisplayName', 'O1' )
semilogy(xref(2:4),ref(2:4,2),'-', 'Color', [152/256 152/256 152/256],'LineWidth',1.0, 'DisplayName', 'O2' )


ylabel('Max Error','fontsize',12)
xlabel('glevel','fontsize',12)
title('Tg Vector Recon Max Errors','fontsize',12 )
box off;

%----------------
%2st graph
%-------------
subplot2=subplot(1,2,2);
set(subplot2, 'Position',[0.52 0.12 width height]);

variable='indmax'

%Reference lines
ref1=0.0025;
ref2=0.04;

hold off;
i=1;
yds=ds(strcmp(ds.grid,grids(i)) , {'Glevel',variable});
x=yds.Glevel(1:end);
y=yds(1:end,2); 
y

semilogy(x,y ,[linestyles{i} Markers(i)],'LineWidth',linewidth, ...
    'MarkerSize', marksize, 'DisplayName', char(names1(i)), 'Color', colors(i,:))
hold on

for i=2:3
yds=ds(strcmp(ds.grid,grids(i)) , {'Glevel',variable});
    x=yds.Glevel(1:end);
    y=yds(1:end,2); 
    y
semilogy(x, y,[linestyles{i} Markers(i)],'LineWidth',linewidth, 'MarkerSize', marksize, 'DisplayName', char(names1(i)), 'Color', colors(i,:))
end
n=size(y);
%Reference lines
ref(1,1)=ref1;
ref(1,2)=ref2;
xref(1)=1;
for j=2:size(y)
    xref(j)=j;
    ref(j, 1)=ref(j-1, 1)/2;
    ref(j, 2)=ref(j-1, 2)/4;
end
semilogy(xref(2:4),ref(2:4,1),'-', 'Color', [152/256 152/256 152/256],'LineWidth',1.0, 'DisplayName', 'O1' )
semilogy(xref(2:4),ref(2:4,2),'-', 'Color', [152/256 152/256 152/256],'LineWidth',1.0, 'DisplayName', 'O2' )


ylabel('Max Error','fontsize',12)
xlabel('glevel','fontsize',12)
title('Max Inconst Index','fontsize',12 )
box off;

leg2=legend(subplot2, 'show'); 
set(leg2, 'Position',[0.84 0.45 0.13 0.13]);
legend boxoff;
