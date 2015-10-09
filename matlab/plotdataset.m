
marcvec=cellstr(char('-*r','-og','-sb','-^c','-+m','-dy','-vk','-pk'));

names=['HCt SCVT avg'; 'HCt HR95 avg'; 'HCm SCVT avg'; 'HCm HR95 avg'; 'HCm SCVT bar'; 'HCm HR95 bar'];
names1=cellstr(names);
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

fig=figure('Color',[1 1 1]);
clf;

set(fig, 'Position', [100 100 800 600])
operator='divuh';
variable='MaxError';
xtitle='Maximum Error';

nmtds=2;
methods=cellstr(char('swm_tc2_HTC_vrectrsk_crectrsk_sinttrsk_grdtrsk', ...
    'swm_tc2_HTC_vrectrsk_crectrsk_sinttrsk_grdtrsk'));

methodnames=cellstr(char('TRSK-HCT-HR95',...
    'TRSK-HCT-SCVT'));

grids=cellstr(char('HR95JT_00',...
    'icos_pol_nopt_'));

subplot1=subplot(2,2,1);
hold off
i=1;
yds=ds(strcmp(ds.Operator,operator) & strcmp(ds.Methods,methods(i)) & ...
   strcmp(ds.Grid,grids(i)) , {'Mesh', 'Glevel',variable});
x=yds.Glevel(2:end);
y=yds.MaxError(2:end);

semilogy(x,y ,[linestyles{1} Markers(1)],'LineWidth',linewidth, ...
    'MarkerSize', marksize, 'DisplayName', char(names1(1)), 'Color', colors(1,:))
hold on
for i=2:nmtds
yds=ds(strcmp(ds.Operator,operator) & strcmp(ds.Methods,methods(i)) & ...
   strcmp(ds.Grid,grids(i)) , {'Mesh', 'Glevel',variable});
    x=yds.Glevel(2:end);
    y=yds.MaxError(2:end);
 semilogy(x, y,[linestyles{i} Markers(i)],'LineWidth',linewidth, 'MarkerSize', marksize, 'DisplayName', char(names1(i)), 'Color', colors(i,:))
end
%Reference lines
ref(1,1)=0.000015;
ref(1,2)=0.00065;
xref(1)=1;
for j=2:size(y)
    xref(j)=j;
    ref(j, 1)=ref(j-1, 1)/2;
    ref(j, 2)=ref(j-1, 2)/4;
end
n=size(y);
semilogy(xref(n-3:n),ref(n-3:n,1),'-', 'Color', [152/256 152/256 152/256],'LineWidth',1.0, 'DisplayName', 'O1' )
semilogy(xref(n-3:n),ref(n-3:n,2),'-', 'Color', [152/256 152/256 152/256],'LineWidth',1.0, 'DisplayName', 'O2' )
ylabel('Max Error','fontsize',12)
xlabel('glevel','fontsize',12)

title('Divergence Maximum Errors','fontsize',12 )
box off;

plot(1)



