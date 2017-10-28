%% Liear analysis
% EEN scheme on an f plane
% P. Peixoto - 2015
%-------------------------
clear;
format shortE 
%close all

%% Parameters        
mtd=0;  % Method used (as delta_B in paper)
        %  0 : original
        %  1 : modified/corrected scheme  
        %  2 : Sadourny Energy conserving scheme 
        %  3 : Sadourny Enstrophy conserving scheme
         

vtc=0;  % Vertical coordinate (as delta_I in paper)
        % 0 for height 
        % 1 for isopycnal 

%% Grid  configuration
% nx, ny : number of x and y grid points
% lx, ly : domain size in x and y
% nx should be preferably power of 2 to speedup code
m=7;
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
par.h0=1;
par.u0=10;
par.v0=0;
par.f0=2/grd.dx;
par.g=1;
%dt=0.00001*grd.dx;

%Non dimensional parameters
par.Ru=2*par.u0/grd.dx/par.f0;    % Rossby Number in u direction 
par.Rv=2*par.v0/grd.dy/par.f0;    % Rossby Number in v direction 
par.c=sqrt(par.g*par.h0);       % Gravity wave speed
par.Fu=par.u0/par.c;            % Froude Number in u direction 
par.Fv=par.v0/par.c;            % Froude Number in v direction 
par.B=par.c^2/par.f0^2/grd.lx^2; %Burgers number
%par.cfl=grd.dt*(u0/grd.dx+v0/grd.dy);      % CFL number
nxy=grd.nx*grd.ny;

var0.u(1:end, 1:end)=par.u0;
var0.v(1:end, 1:end)=par.v0;
var0.h(1:end, 1:end)=par.h0;

var = var0;
% var1=timestep(0, pert, var0, grd, par, 0);
[tu0, tv0, th0, var]=tendencies(var, grd, par);


kappa(grd.nx/2+1)=0;
lambda(grd.ny/2+1)=0;
gr(grd.nx/2+1, grd.ny/2+1)=0;

%i=1;
%ikappa=2;
for i=1:grd.nx/2+1
    ikappa=i-1;
    kappa(i)=2*pi*ikappa/grd.nx;
    disp(kappa(i))
    for j=1:grd.ny/2+1
        jlambda=j-1;
        lambda(j)=2*pi*jlambda/grd.ny;
        gr(i, j)=growth_rate_linanalysis(ikappa, jlambda, grd, par, var0, tu0, tv0, th0);
    end
end


fig=figure('Color',[1 1 1]);
clf;
contourf(kappa, lambda, gr', 100, 'LineColor','none')
xlabel('\kappa')
ylabel('\lambda')
colorbar
caxis([0,3])
