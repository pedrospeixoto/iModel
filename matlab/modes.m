%% Calculate normal modes of SWM 
% matrix read from files in sparse form

readdata=0;
if readdata==1
    tic;    
        clear;
    %Numerical modes
     B=dlmread('swm_tc11_HTC_vrectrsk_crectrsk_sinttrsk_grdtrsk_linearmatrixsparse_HR95JT_004.txt');
     A=sparse(B(:,1),B(:,2), B(:,3));
     e11trsk=eigs(A, 10240);   %%sizeofmatriz -2 to be able to use sparse matrix - matlab limitation
     1
     toc
     tic;
     B=dlmread('swm_tc12_HTC_vrectrsk_crectrsk_sinttrsk_grdtrsk_linearmatrixsparse_HR95JT_004.txt');
     A=sparse(B(:,1),B(:,2), B(:,3));
     e12trsk=eigs(A, 10240);   %%sizeofmatriz -2 to be able to use sparse matrix - matlab limitation
     2
     toc
     tic;
     B=dlmread('swm_tc11_HC_vrecperhx_crecpered_sintbary_grdtrsk_linearmatrixsparse_HR95JT_004.txt');
     A=sparse(B(:,1),B(:,2), B(:,3));
     e11new=eigs(A, 10240);   %%sizeofmatriz -2 to be able to use sparse matrix - matlab limitation
     3
     toc
     tic;
     B=dlmread('swm_tc12_HC_vrecperhx_crecpered_sintbary_grdtrsk_linearmatrixsparse_HR95JT_004.txt');
     A=sparse(B(:,1),B(:,2), B(:,3));
     e12new=eigs(A, 10240);   %%sizeofmatriz -2 to be able to use sparse matrix - matlab limitation
     4
     toc
     tic;
     clearvars A B
end

%Analytical modes
ndof=10242;
ncells=2562;
nvert_tr=5120;
nedges=7680;
ngeo=5120;
ningrav=2*ncells;
gha2=2.46351265002725E-009;
f=1.46E-004;
ngeo2=(ngeo-2)/2; %two zero eigenvaleus were not calculated - so we will take them of the geo modes
ningrav2=ningrav/2;
nmodes2=ngeo2+ningrav2;
first=10240/2;

%Analytical modes
k=0;
for i=0:50
    for j=1:i*2+1
        k=k+1;
        ingrav_exact(k)=sqrt(f^2+i*(i+1)*gha2);
    end
end

e11trsk_s=sort(imag(e11trsk));
e12trsk_s=sort(imag(e12trsk));
e11new_s=sort(imag(e11new));
e12new_s=sort(imag(e12new));

geo_trsk11   =e11trsk_s(first:first+ngeo2-2);
ingrav_trsk11=e11trsk_s(first+ngeo2:end);
geo_trsk12   =e12trsk_s(first:first+ngeo2-2);
ingrav_trsk12=e12trsk_s(first+ngeo2:end);
geo_new11    =e11new_s (first:first+ngeo2-2);
ingrav_new11 =e11new_s (first+ngeo2:end);
geo_new12    =e12new_s(first:first+ngeo2-2);
ingrav_new12 =e12new_s(first+ngeo2:end);

names1=cellstr(char('exact',...
    'trsk-f', ...
    'trsk-r', ...
    'modf-f', ...
    'modf-r'));


fig=figure('Color',[1 1 1]);
clf;
set(fig, 'Position', [100 100 800 600])

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
plot(ingrav_exact ,[linestyles{1} Markers(1)],'LineWidth',linewidth, ...
    'MarkerSize', marksize, 'DisplayName', char(names1(1)), 'Color', colors(1,:))
hold on
plot(ingrav_trsk11 ,[linestyles{2} Markers(1)],'LineWidth',linewidth, ...
    'MarkerSize', marksize, 'DisplayName', char(names1(2)), 'Color', colors(2,:))
plot(ingrav_trsk12 ,[linestyles{3} Markers(1)],'LineWidth',linewidth, ...
    'MarkerSize', marksize, 'DisplayName', char(names1(3)), 'Color', colors(3,:))
plot(ingrav_new11 ,[linestyles{4} Markers(1)],'LineWidth',linewidth, ...
    'MarkerSize', marksize, 'DisplayName', char(names1(4)), 'Color', colors(4,:))
plot(ingrav_new12 ,[linestyles{5} Markers(1)],'LineWidth',linewidth, ...
    'MarkerSize', marksize, 'DisplayName', char(names1(5)), 'Color', colors(5,:))
xlabel('','fontsize',12)
ylabel('Frequency','fontsize',12)
title('Inertial-Gravity Modes','fontsize',12 )
box off;
leg1=legend(subplot1, 'show'); 
set(leg1, 'Position',[0.85 0.65 0.13 0.13],'fontsize',12);
legend boxoff;
hold off
hold off

subplot2=subplot(2,1,2);
set(subplot2, 'Position',[0.08 0.08 0.80 0.38]);
hold off;
plot(geo_trsk11 ,[linestyles{2} Markers(1)],'LineWidth',linewidth, ...
    'MarkerSize', marksize, 'DisplayName', char(names1(2)), 'Color', colors(2,:))
hold on
plot(geo_trsk12 ,[linestyles{3} Markers(1)],'LineWidth',linewidth, ...
    'MarkerSize', marksize, 'DisplayName', char(names1(3)), 'Color', colors(3,:))
plot(geo_new11 ,[linestyles{4} Markers(1)],'LineWidth',linewidth, ...
    'MarkerSize', marksize, 'DisplayName', char(names1(4)), 'Color', colors(4,:))
plot(geo_new12 ,[linestyles{5} Markers(1)],'LineWidth',linewidth, ...
    'MarkerSize', marksize, 'DisplayName', char(names1(5)), 'Color', colors(5,:))
xlabel('Modes','fontsize',12)
ylabel('Frequency','fontsize',12)
title('Geostrophic Modes','fontsize',12 )
box off;
leg2=legend(subplot2, 'show'); 
set(leg2, 'Position',[0.85 0.25 0.13 0.13],'fontsize',12);
legend boxoff;
hold off


