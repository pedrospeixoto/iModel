%% Calculate normal modes of SWM 
% matrix read from files in sparse form

readdata=0;
if readdata==1
    tic;    
    clear;
    %Numerical modes
     B=dlmread('swm_tc32_HTC_vrectrsk_crectrsk_sinttrsk_grdtrsk_hol1.0_linearmatrixsparse_HR95JT_003.txt');
     A=sparse(B(:,1),B(:,2), B(:,3));
    e1=eigs(A, 2560);   %%sizeofmatriz -2 to be able to use sparse matrix - matlab limitation
    
     B=dlmread('swm_tc32_HTC_vrectrsk_crectrsk_sinttrsk_grdtrsk_hol10.0_linearmatrixsparse_HR95JT_003.txt');
     A=sparse(B(:,1),B(:,2), B(:,3));
     e10=eigs(A, 2560);   %%sizeofmatriz -2 to be able to use sparse matrix - matlab limitation
     
     clearvars A B
end

%Analytical modes
ndof=2562;
ncells=642;
nvert_tr=1280;
nedges=1920;
ngeo=5120;
ningrav=2*ncells;
gha2=2.46351265002725E-009;
f=1.46E-004;
ngeo2=(ngeo-2)/2; %two zero eigenvaleus were not calculated - so we will take them of the geo modes
ningrav2=ningrav/2;
nmodes2=ngeo2+ningrav2;
first=10240/2;


[e1_rmax, e1_rmaxi]=max(real(e1))
e1_i=imag(e1);
e1_compfreq=e1_i(e1_rmaxi)


e1_i=sort(imag(e1));
e1_r=sort(real(e1));
e10_i=sort(imag(e10));
e10_r=sort(real(e10));

fig=figure('Color',[1 1 1]);
clf;
set(fig, 'Position', [100 100 800 600])

names1=cellstr(char('h1-real',...
    'h10-real', ...
    'h10-imag', ...
    'h10-real' ...
    ));

marcvec=cellstr(char('-*r','-og','-sb','-^c','-+m','-dy','-vk','-pk'));

marksize=6;
linewidth=1.8;
linestyles = cellstr(char('-','--','-.','--','-.','--','-.',':','-',':','-',':',...
    '-.','--','-',':','-.','--','-',':','-.'));

Markers=[' ', 'o','s','+','x','^','d','v','^','<','>','p','h','.',...
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


subplot1=subplot(2,1,1);
set(subplot1, 'Position',[0.08 0.54 0.80 0.40]);
hold off;
plot(e1_r ,[linestyles{1} Markers(1)],'LineWidth',linewidth, ...
    'MarkerSize', marksize, 'DisplayName', char(names1(1)), 'Color', colors(1,:))
hold on
plot(e10_r ,[linestyles{2} Markers(1)],'LineWidth',linewidth, ...
    'MarkerSize', marksize, 'DisplayName', char(names1(2)), 'Color', colors(2,:))
xlabel('','fontsize',12)
ylabel('Frequency','fontsize',12)
title('Modes depth=1','fontsize',12 )
box off;
leg1=legend(subplot1, 'show'); 
set(leg1, 'Position',[0.85 0.65 0.13 0.13],'fontsize',12);
legend boxoff;
hold off
hold off


subplot2=subplot(2,1,2);
set(subplot2, 'Position',[0.08 0.08 0.80 0.38]);
hold off;
plot(e1_i ,[linestyles{2} Markers(1)],'LineWidth',linewidth, ...
    'MarkerSize', marksize, 'DisplayName', char(names1(3)), 'Color', colors(1,:))
hold on
plot(e10_i ,[linestyles{3} Markers(1)],'LineWidth',linewidth, ...
    'MarkerSize', marksize, 'DisplayName', char(names1(4)), 'Color', colors(2,:))
xlabel('Modes','fontsize',12)
ylabel('Frequency','fontsize',12)
title('Modes depth=10','fontsize',12 )
box off;
leg2=legend(subplot2, 'show'); 
set(leg2, 'Position',[0.85 0.25 0.13 0.13],'fontsize',12);
legend boxoff;
hold off


fig=figure('Color',[1 1 1]);
clf;
set(fig, 'Position', [100 100 800 600])

names1=cellstr(char('h1',...
    'h10' ...
    ));

subplot1=subplot(2,1,1);
set(subplot1, 'Position',[0.08 0.54 0.80 0.36]);

hold off;
plot(e1 , 'o')
hold on
xlabel('Real','fontsize',12)
ylabel('Imag','fontsize',12)
title('Modes depth=1','fontsize',12 )
box off;
hold off

subplot2=subplot(2,1,2);
set(subplot2, 'Position',[0.08 0.08 0.80 0.36]);
hold off;
plot(e10 , 'o')
hold on
xlabel('Real','fontsize',12)
ylabel('Imag','fontsize',12)
title('Modes depth=10','fontsize',12 )
box off;
hold off



