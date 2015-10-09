%% read_data %call read_data.m

% import data from tc2.txt as a dataset named ds

fig=figure('Color',[1 1 1]);
clf;

set(fig, 'Position', [100 100 800 300])

%General parameters
nmtds=4;
methods=cellstr(char('swm_tc2_nt10000_HTC_vrectrsk_crectrsk_sinttrsk_grdtrsk', ...
    'swm_tc2_nt10000_HTC_vrectrsk_crectrsk_sinttrsk_grdtrsk', ...
    'swm_tc2_nt10000_HC_vrecperhx_crecpered_sintbary_grdtrsk', ...
    'swm_tc2_nt10000_HC_vrecperhx_crecpered_sintbary_grdtrsk'));

methodnames=cellstr(char('TRSK-HCT-SCVT',...
    'TRSK-HCT-HR95', ...
    'MODF-HCM-SCVT', ...
    'MODF-HCM-HR95'));

grids=cellstr(char('icos_pol_scvt_h1_',...
    'HR95JT_00', ...
    'icos_pol_scvt_h1_',...
    'HR95JT_00'));

height=0.70;
width=0.40;

%----------------
%1st graph
%-------------
subplot1=subplot(1,2,1);
set(subplot1, 'Position',[0.08 0.2 width height]);

time=12
variable='errormax_h'

%Reference lines
ref1=0.0001;
ref2=0.01;
plotgraph_tc2(ds, time, variable, methods, nmtds, grids, methodnames, ref1, ref2)
ylabel('Max Error','fontsize',12)
xlabel('glevel','fontsize',12)
title('TC2 Max h Errors','fontsize',12 )
box off;
axis([3 9 0.0000001 0.01])

%2nd graph
subplot2 = subplot(1,2,2);
set(subplot2, 'Position',[0.57 0.2 width height]);
variable='error2_h'
%set(subplot2, 'Position',[0.48 0.58 0.30 0.35]);
ref1=0.00009;
ref2=0.006;
plotgraph_tc2(ds, time, variable, methods, nmtds, grids, methodnames, ref1, ref2)
ylabel('RMS Error','fontsize',12)
xlabel('glevel','fontsize',12)
title('TC2 RMS Errors','fontsize',12 )
box off;
%leg2=legend(subplot2, 'show'); 
%set(leg2, 'Position',[0.82 0.75 0.13 0.13]);
%legend boxoff;
hold off
axis([3 9 0.0000001 0.001])

leg=legend(subplot2, 'show'); 
set(leg, 'Position',[0.09 0.02 0.8 0.03],'Orientation','horizontal');
legend boxoff;
hold off


