%% read_data %call read_data.m

fig=figure('Color',[1 1 1]);
clf;

set(fig, 'Position', [100 100 800 300])

%General parameters
nmtds=4;
methods=cellstr(char('swm_tc32_dt50_HTC_vrectrsk_crectrsk_sinttrsk_grdtrsk_hol100.0', ...
    'swm_tc32_dt50_HTC_vrectrsk_crectrsk_sinttrsk_grdtrsk_hol100.0', ...
    'swm_tc32_dt50_HC_vrecperhx_crecpered_sintbary_grdtrsk_hol100.0', ...
    'swm_tc32_dt50_HC_vrecperhx_crecpered_sintbary_grdtrsk_hol100.0' ...
    ));

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
ref1=0.003;
ref2=0.1;
plotgraph_tc32(ds, time, variable, methods, nmtds, grids, methodnames, ref1, ref2)
ylabel('Max Error','fontsize',12)
xlabel('glevel','fontsize',12)
title('Thin layer Max Errors (h)','fontsize',12 )
box off;

%2nd graph
subplot2 = subplot(1,2,2);
set(subplot2, 'Position',[0.57 0.2 width height]);
variable='error2_h'
%set(subplot2, 'Position',[0.48 0.58 0.30 0.35]);
ref1=0.001;
ref2=0.06;
plotgraph_tc32(ds, time, variable, methods, nmtds, grids, methodnames, ref1, ref2)
ylabel('RMS Error','fontsize',12)
xlabel('glevel','fontsize',12)
title('Thin layer RMS Errors (h)','fontsize',12 )
box off;
%leg2=legend(subplot2, 'show'); 
%set(leg2, 'Position',[0.82 0.75 0.13 0.13]);
%legend boxoff;
hold off

leg=legend(subplot2, 'show'); 
set(leg, 'Position',[0.09 0.02 0.8 0.03],'Orientation','horizontal');
legend boxoff;
hold off


