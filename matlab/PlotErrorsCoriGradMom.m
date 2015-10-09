%% read_data %call read_data.m

fig=figure('Color',[1 1 1]);
clf;

set(fig, 'Position', [100 100 800 800])

%General parameters
nmtds=4;
methods=cellstr(char('swm_tc2_HTC_vrectrsk_crectrsk_sinttrsk_grdtrsk', ...
    'swm_tc2_HTC_vrectrsk_crectrsk_sinttrsk_grdtrsk',...
    'swm_tc2_HC_vrecperhx_crecpered_sintbary_grdtrsk',...
    'swm_tc2_HC_vrecperhx_crecpered_sintbary_grdtrsk'));

methodnames=cellstr(char('TRSK-HCT-SCVT',...
    'TRSK-HCT-HR95', ...
    'MODF-HCM-SCVT', ...
    'MODF-HCM-HR95'));

grids=cellstr(char('icos_pol_scvt_h1_',...
    'HR95JT_00', ...
    'icos_pol_scvt_h1_',...
    'HR95JT_00'));

height=0.25;
width=0.34;

%----------------
%1st graph
%-------------
subplot1=subplot(3,2,1);
set(subplot1, 'Position',[0.08 0.72 width height]);

operator='uhq_perp'
variable='MaxError'

%Reference lines
ref1=0.000006;
ref2=0.00025;
plotgraph(ds, operator, variable, methods, nmtds, grids, methodnames, ref1, ref2)
ylabel('Max Error','fontsize',12)
%xlabel('glevel','fontsize',12)
title('Perp Max Errors','fontsize',12 )
box off;

%2nd graph
subplot2 = subplot(3,2,2);
set(subplot2, 'Position',[0.55 0.72 width height]);
operator='uhq_perp'
variable='RMSError'
%set(subplot2, 'Position',[0.48 0.58 0.30 0.35]);
ref1=0.0009;
ref2=0.06;
plotgraph(ds, operator, variable, methods, nmtds, grids, methodnames, ref1, ref2)
ylabel('RMS Error','fontsize',12)
%xlabel('glevel','fontsize',12)
title('Perp RMS Errors','fontsize',12 )
box off;
%leg2=legend(subplot2, 'show'); 
%set(leg2, 'Position',[0.82 0.75 0.13 0.13]);
%legend boxoff;
hold off


%----------------
%3rd graph
%-------------
subplot3=subplot(3,2,3);
set(subplot3, 'Position',[0.08 0.40 width height]);

operator='grad_ghbK'
variable='MaxError'

%Reference lines
ref1=0.0000015;
ref2=0.000095;
plotgraph(ds, operator, variable, methods, nmtds, grids, methodnames, ref1, ref2)
ylabel('Max Error','fontsize',12)
%xlabel('glevel','fontsize',12)
title('Grad Max Errors','fontsize',12 )
box off;

%4th graph
subplot4 = subplot(3,2,4);
set(subplot4, 'Position',[0.55 0.4 width height]);
operator='grad_ghbK'
variable='RMSError'
ref1=0.00009;
ref2=0.006;
plotgraph(ds, operator, variable, methods, nmtds, grids, methodnames, ref1, ref2)
ylabel('RMS Error','fontsize',12)
%xlabel('glevel','fontsize',12)
title('Grad RMS Errors','fontsize',12 )
box off;
%leg4=legend(subplot4, 'show'); 
%set(leg4, 'Position',[0.82 0.45 0.13 0.13]);
%legend boxoff;
hold off

%5th graph
subplot5 = subplot(3,2,5);
set(subplot5, 'Position',[0.08 0.08 width height]);
operator='momeq'
variable='MaxError'
ref1=0.000003;
ref2=0.00015;
plotgraph(ds, operator, variable, methods, nmtds, grids, methodnames, ref1, ref2)
ylabel('Max Error','fontsize',12)
xlabel('glevel','fontsize',12)
title('Mom Eq Max Errors','fontsize',12 )
box off;

%6th graph
subplot6 = subplot(3,2,6);
set(subplot6, 'Position',[0.55 0.08  width height]);
operator='momeq'
variable='RMSError'
ref1=0.000001;
ref2=0.00005;
plotgraph(ds, operator, variable, methods, nmtds, grids, methodnames, ref1, ref2)
ylabel('RMS Error','fontsize',12)
xlabel('glevel','fontsize',12)
title('Mom Eq RMS Errors','fontsize',12 )
box off;
leg6=legend(subplot6, 'show'); 
set(leg6, 'Position',[0.09 0.00 0.8 0.03],'Orientation','horizontal');
legend boxoff;
hold off


