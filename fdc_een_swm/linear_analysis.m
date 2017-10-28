%% Liear analysis
% EEN scheme on an f plane
% P. Peixoto - 2015
%-------------------------
clear;
format shortE 
close all

%% Parameters        
mtd=0;  % Method used (as delta_B in paper)
        %  0 : original
        %  1 : modified/corrected scheme  
        %  2 : Sadourny Energy conserving scheme 
        %  3 : Sadourny Enstrophy conserving scheme
         

vtc=1;  % Vertical coordinate (as delta_I in paper)
        % 0 for height 
        % 1 for isopycnal 

%% Grid  configuration
% nx, ny : number of x and y grid points
% lx, ly : domain size in x and y
% nx should be preferably power of 2 to speedup code
m=4;
nx=2^m;
ny=2^m;
lx=1;
ly=1;
 
tic

%% Initialize grid
grd=initialize_grid(nx, ny, lx, ly, 1, 1);


%% Initialize test case - just to get structure right
[var0, par, varmax]=initialize_tc(grd, 1, mtd, vtc);

%Overwrite initial conditions
par.h0=2.5/16^2;
par.u0=0.5*100/16;
par.v0=10;
par.f0=10;
par.g=10;
pert=0.00000001;
%dt=0.00001*grd.dx;
%Non dimansional parameters
par.Ru=2*par.u0/grd.dx/par.f0;    % Rossby Number in u direction 
par.Rv=2*par.v0/grd.dy/par.f0;    % Rossby Number in v direction 
par.c=sqrt(par.g*par.h0);       % Gravity wave speed
par.Fu=par.u0/par.c;            % Froude Number in u direction 
par.Fv=par.v0/par.c;            % Froude Number in v direction 
par.B=par.c^2/par.f0^2/grd.lx^2; %Burgers number
%par.cfl=grd.dt*(u0/grd.dx+v0/grd.dy);      % CFL number

var0.u(1:end, 1:end)=par.u0;
var0.v(1:end, 1:end)=par.v0;
var0.h(1:end, 1:end)=par.h0;

var = var0;
% var1=timestep(0, pert, var0, grd, par, 0);
[tu0, tv0, th0, var]=tendencies(var, grd, par);


nxy=grd.nx*grd.ny;
h0=reshape(var0.h, nxy, 1);
u0=reshape(var0.u, nxy, 1);
v0=reshape(var0.v, nxy, 1);

M(3*nxy, 3*nxy)=0.;

k=0;
%Perturb h
for j=1:grd.ny
    for i=1:grd.nx
        k=k+1;
        var=var0;
        var.h(i, j)=var0.h(i,j)+pert;
        %var=timestep(0, dt, var, grd, par, 0);
        %h=(var.h-var1.h)/pert;
        %u=(var.u-var1.u)/pert;
        %v=(var.v-var1.v)/pert;
        % And save as column of system matrix for h, u and v
        %M(1:nxy        ,k) = reshape(h, nxy, 1);
        %M(nxy+1:2*nxy  ,k) = reshape(u, nxy, 1);
        %M(2*nxy+1:3*nxy,k) = reshape(v, nxy, 1);
        [tu, tv, th, var]=tendencies(var, grd, par);
        M(1:nxy        ,k) = reshape(th-th0, nxy, 1)/pert;
        M(nxy+1:2*nxy  ,k) = reshape(tu-tu0, nxy, 1)/pert;
        M(2*nxy+1:3*nxy,k) = reshape(tv-tv0, nxy, 1)/pert;
    end
end


%Perturb u
for j=1:grd.ny
    for i=1:grd.nx
        k=k+1;
        var=var0;
        var.u(i, j)=var0.u(i,j)+pert;
        %var=timestep(0, dt, var, grd, par, 0);
        %h=(var.h-var1.h)/pert;
        %u=(var.u-var1.u)/pert;
        %v=(var.v-var1.v)/pert;
        % And save as column of system matrix for h, u and v
        %M(1:nxy        ,k) = reshape(h, nxy, 1);
        %M(nxy+1:2*nxy  ,k) = reshape(u, nxy, 1);
        %M(2*nxy+1:3*nxy,k) = reshape(v, nxy, 1);
        [tu, tv, th, var]=tendencies(var, grd, par);
        M(1:nxy        ,k) = reshape(th-th0, nxy, 1)/pert;
        M(nxy+1:2*nxy  ,k) = reshape(tu-tu0, nxy, 1)/pert;
        M(2*nxy+1:3*nxy,k) = reshape(tv-tv0, nxy, 1)/pert;
    end
end

%Perturb v
for j=1:grd.ny
    for i=1:grd.nx
        k=k+1;
        var=var0;
        var.v(i, j)=var0.v(i,j)+pert;
        %var=timestep(0, dt, var, grd, par, 0);
        %h=(var.h-var1.h)/pert;
        %u=(var.u-var1.u)/pert;
        %v=(var.v-var1.v)/pert;        
        % And save as column of system matrix for h, u and v
        %M(1:nxy        ,k) = reshape(h, nxy, 1);
        %M(nxy+1:2*nxy  ,k) = reshape(u, nxy, 1);
        %M(2*nxy+1:3*nxy,k) = reshape(v, nxy, 1);
        [tu, tv, th, var]=tendencies(var, grd, par);
        M(1:nxy        ,k) = reshape(th-th0, nxy, 1)/pert;
        M(nxy+1:2*nxy  ,k) = reshape(tu-tu0, nxy, 1)/pert;
        M(2*nxy+1:3*nxy,k) = reshape(tv-tv0, nxy, 1)/pert;
    end
end
toc

% %Analytical inertia-gravity frequencies for u0=0
% ingrav(1:grd.nx, 1:grd.ny)=0.
% for j=1:grd.ny
%     for i=1:grd.nx
%         ingrav(i,j)=sqrt(par.f0^2+par.c^2*(j^2+i^2));
%     end 
% end
% 
% plot(ingrav)
% error

%A=sparse(M);
tic;
%z=eig(M);
%s=sign(z).*log(z)/pert;
s = eig(M);
%s=eigs(A, nxy-1);
[si,ix]=sort(imag(s));
%si=sort(imag(s));
sr = real(s(ix));

%Adjust the notation with the paper (divide by f0)
sr=sr/par.f0;
si=si/par.f0;

toc


fig=figure('Color',[1 1 1], 'Position', [100, 100, 600, 400]);

width=0.8;
height=0.3;
sub1=subplot(2,1,1);
set(sub1, 'Position',[0.1 0.65 width height]);

plot(-3*nxy/2:1:3*nxy/2-1, si);
title('Imag part of eigenvalues (stable)')
ylabel('Frequency')
%xlabel('Mode index')

sub2=subplot(2,1,2);
set(sub2, 'Position',[0.1 0.2 width height]);
plot(-3*nxy/2:1:3*nxy/2-1, sr);
title('Real part of eigenvalues (unstable)')
ylabel('Growth rate')
xlabel('Mode index')

coment=annotation('textbox',[0.05 0.01 1 0.05]);
set(coment,'linestyle','none')
set(coment,'FitBoxToText','on')

set(coment,'string',[...
    ' R_u= ',num2str(par.Ru), ...
    '     F_u= ',num2str(par.Fu), ...
    '     Mtd = ', num2str(par.mtd), ...
    '     Vtc = ', num2str(par.vtc), ...
    '     N =', num2str(grd.nx), ...
    ], 'fontsize',11);
