%% Shallow Water Model 2d
% EEN scheme on an f plane
% P. Peixoto - 2015
%-------------------------
%clear;
format shortE 
close all

%% Parameters
tc=5;   % Test case:
        %   0 : Read data from file
        %   1 : u=sin(2pi y), h to balance
        %   2 : v=sin(2pi x), h to balance
        %   3 : Gaussian 'dam break'
        %   4 : Same as 1, but balanced state put in 'b' and h small constant
        %   5 : Hollingsworth test case - similar to 4, but with specially chosen parameters
        %   6 : Constant zonal velocity with forcing term to balance
        %   7 : Advection test case - for validation of timestepping scheme (not debugged)
        %   8 : Test case 5 set up to catch nonlinear instability
        %   9 : Steady state similar to 1 but to match sweet notation
        
mtd=0;  % Method used (as delta_B in paper)
        %  0 : original
        %  1 : modified/corrected scheme  
        %  2 : Sadourny Energy conserving scheme 
        %  3 : Sadourny Enstrophy conserving scheme
        
vtc=0;  % Vertical coordinate (as delta_I in paper)
        % 0 for height (uses abs vorticity in rot term)
        % 1 for isopycnal (used PV in rot term)

%% Grid  configuration
% nx, ny : number of x and y grid points
% lx, ly : domain size in x and y
% nx should be preferably power of 2 to speedup code
m=7;
nx=2^m; %Could be small if test case only depends on y
ny=2^m;
lx=1;
ly=1;
 

%% Time stepping parameters
Tmax=1;        % Final time
dt=0.1*1/1000;    %Tmax/nt;    % Time step
nt=Tmax/dt;         % Number of time steps
Tstop=1;    % Intermediate time to stop
graphtime=10*Tmax/(nt);  % Time between graph plots

tic

%% Initialize grid
grd=initialize_grid(nx, ny, lx, ly, nt, Tmax);


%% Initialize test case
[var0, par, varmax]=initialize_tc(grd, tc, mtd, vtc);

%Save initial conditions
var=var0;

%Calculate tendencies just to get the derived quantities in var
[tu, tv, th, var]=tendencies(var, grd, par);

%plot_var(var.gy+var.zu, grd, par, 'h', 'h');


%% Loop over time
%Time stepping
t=0;

fig=figure('Color',[1 1 1], 'Position', [100, 100, 800, 800]);
%fig1=figure('Color',[1 1 1], 'Position', [900, 100, 600, 800]);
fig2=figure('Color',[1 1 1], 'Position', [900, 600, 600, 300]);
fig3=figure('Color',[1 1 1], 'Position', [900, 100, 600, 300]);

%plot_vars_t(t, var, grd, par, fig);
%plot_var_t(t, var.u, grd, par, 'u', 'u', fig);
%pause()

kenergy0=calc_energy(var, grd, par);
enstrophy0=calc_enstrophy(var, grd, par);

energy(1:nt+1)=0.;
var.gx=calc_gradx(var.ke, grd);
var.gy=calc_grady(var.ke, grd);
max_abs_gradke(1:nt+1)=0.;
max_abs_gradke(1)=max(max(sqrt(var.gx.*var.gx+var.gy.*var.gy)));

energy(1)=kenergy0;
enstrophy(1:nt+1)=0.;
herror(1:nt)=0.;

%profile on
for k=1:nt %445 
    t=dt*k;
    
    var=timestep(t, dt, var, grd, par, 1);

    %energy(k)=(calc_energy(var, grd, par)-kenergy0)/kenergy0;
    energy(k+1)=calc_energy(var, grd, par); %-kenergy0)/kenergy0;
    enstrophy(k+1)=(calc_enstrophy(var, grd, par)-enstrophy0)/enstrophy0;
    var.gy=abs(calc_grady(var.v, grd));
    max_abs_gradke(k+1)=max(max(var.gy));
    
    if(mod(t,graphtime)==0)
        t
        %plot_vars_t(t, var, grd, par, fig);
        %plot_var_t(t, var.h-var0.h, grd, par, 'h', 'h error ', fig2);
        %plot_var_t(t, var.h, grd, par, 'h', 'h', fig);
        %varmax=plot_yslice(t, k, var, var0, varmax, grd, fig);   
        varmax=plot_yslice4(t, k, var, var0, varmax, grd, fig);   
        plot_en_evol(t, k, energy, max_abs_gradke, fig2);   
        herror=plot_evol(t, k, var, var0, herror, grd, fig3);   
        %pause();
    end
    if t>Tstop  || herror(k)>100
        break
    end
end 
%profile viewer

%% Outputs
error_max_h=max(max(var.h-var0.h))
error_2_h=sum(sum((var.h-var0.h).*(var.h-var0.h)))/(grd.nx*grd.ny)
%plot_var(var.h-var0.h, grd, par, 'h', 'error');

toc

%Usefull for growth calculation
% log(herror(453)/herror(400))/(par.f0*(time(453)-time(400)))

%Useful for error comparisons
% > semilogy(time(1:n-1080), herror_c0v0(1:n-1080), time(1:n), herror_c0v1(1:n), time(1:n), herror_c1v0(1:n), time(1:n), herror_c1v1(1:n))
% >> title('SWM test case')
% >> axis([0 1.5 0.0000001 1])
% >> axis([0 1.5 0.00001 1])
% >> axis([0 1.5 0.0001 1])
% >> axis([0 1.5 0.0001 0.1])
% >> xlabel('Time')
% >> ylabel('Max Error')
% >> legend('Orig-Height', 'Orig-Isop', 'Modf-Height', 'Modf-Isop')