function plotgraph(ds, operator, variable, methods, nmtds, grids, names1, ref1, ref2)
format short e
marcvec=cellstr(char('-*r','-og','-sb','-^m','-+c','-dy','-vk','-pk'));

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
colors(3,:)=[150/256 150/256 0/256];
colors(4,:)=[0/256 0/256 153/256];
colors(5,:)=[162/256 0/256 162/256];
colors(6,:)=[162/256 162/256 0/256];
colors(7,:)=[192/256 192/256 192/256];

hold off;
i=1;
yds=ds(strcmp(ds.Operator,operator) & strcmp(ds.Methods,methods(i)) & ...
   strcmp(ds.Grid,grids(i)) , {'Mesh', 'Glevel',variable});
x=yds.Glevel(4:end);
y=yds(4:end,3); 
y

semilogy(x,y ,[linestyles{i} Markers(i)],'LineWidth',linewidth, ...
    'MarkerSize', marksize, 'DisplayName', char(names1(i)), 'Color', colors(i,:))
hold on
for i=2:nmtds
yds=ds(strcmp(ds.Operator,operator) & strcmp(ds.Methods,methods(i)) & ...
   strcmp(ds.Grid,grids(i)) , {'Mesh', 'Glevel',variable});
    x=yds.Glevel(4:end);
    y=yds(4:end,3); 
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
semilogy(xref(5:6),ref(5:6,1),'-', 'Color', [152/256 152/256 152/256],'LineWidth',1.0, 'DisplayName', 'O1' )
semilogy(xref(5:6),ref(5:6,2),'-', 'Color', [152/256 152/256 152/256],'LineWidth',1.0, 'DisplayName', 'O2' )

end