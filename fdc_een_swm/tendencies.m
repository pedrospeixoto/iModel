%% Calculate the RHS of the SW equations
% Receives as input u, v, h inside var structure
function [tu, tv, th, var]=tendencies(var, grd, par)

%Calculate Vorticity
var.z=calc_vorticity(var.u, var.v, grd);
var.z=var.z+par.f0;

%Reconstruct u and v at h points
%[uh, vh]=calc_uv_h(u, v, nx, ny);

%Calculate h at u, v and q points
[h_u, h_v, h_q]=calc_h_uvq(var.h, grd, par.mtd);

%Calculate pv and fluxes
var.q=var.z./h_q;
us=var.u.*h_u;
vs=var.v.*h_v;

% Calculate Coriolis term depending on the vertical coordinate
if par.vtc == 1 % Isopycnal coordinates
        
    % Interpolate potential vorticity to velocity points
    [q_y, q_x]=calc_q_uv(var.q, grd);
    
    
    %Calculate Coriolis terms
    var.zv=calc_coriolis_u(var.q, q_y, q_x, vs, grd, par.mtd);
    var.zu=calc_coriolis_v(var.q, q_y, q_x, us, grd, par.mtd);
    
elseif par.vtc == 0  %Height coordinates
    
    % Interpolate vorticity to velocity points
    [z_y, z_x]=calc_q_uv(var.z, grd);
    
    %Calculate Coriolis terms
    var.zv=calc_coriolis_u(var.z, z_y, z_x, var.v, grd, par.mtd);
    var.zu=calc_coriolis_v(var.z, z_y, z_x, var.u, grd, par.mtd);
    
end

%Calculate the Kinetic Energy
var.ke=calc_ke(var.u, var.v, grd, par.mtd);

%Check for neutral non advection term (need to take f0=0)
var.gx=calc_gradx(var.ke, grd);
var.gy=calc_grady(var.ke, grd);
var.tmpx=-(var.zv)+var.gx;
var.tmpy=(var.zu)+var.gy;

%Calculate the bernoulli potential
phi=par.g*(var.h+var.b)+var.ke;

%Calculate gradient terms
var.gx=calc_gradx(phi, grd);
var.gy=calc_grady(phi, grd);

% %Calculate h at u, v and q points
% [h_u, h_v, h_q]=calc_h_uvq(var.h, grd, 0);
%     
% %Calculate pv and fluxes
% %var.q=var.z./h_q;
% us=var.u.*h_u;
% vs=var.v.*h_v;
    
%Calculate divergence term
var.dv=calc_div(us, vs, grd);

%Tendencies
%var.u
tu=var.zv-var.gx;
%var.v
tv=-var.zu-var.gy;
%var.h
th=-var.dv;

%add forcing term for test case 6
if par.tc==6
    tv=tv+par.f0*par.u0;
    var.tmpy=(var.zu)+var.gy-par.f0*par.u0;
end
%pause()
%Returns var with diagnosed quatities as well

end