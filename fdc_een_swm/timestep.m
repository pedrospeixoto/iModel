%% Runge- Kutta 44 time step
% Time stepping routine 
% t : time for which the step is to advanced to (from t-dt -> t)
function varf=timestep(t, dt, var, grd, par, rk4)

vartmp=var;
varf=var;

if rk4==0 % use Euler method
    
    t0=t-dt;
    [tu0, tv0, th0, vartmp]=tendencies(vartmp, grd, par);
    
    t1 = t0 + dt/2;
    varf.u = var.u + dt * tu0;
    varf.v = var.v + dt * tv0;
    varf.h = var.h + dt * th0;
    
    
elseif rk4>0 %use rk44
    
    t0=t-dt;
    [tu0, tv0, th0, vartmp]=tendencies(vartmp, grd, par);
    
    t1 = t0 + dt/2;
    vartmp.u = var.u + dt * tu0/2;
    vartmp.v = var.v + dt * tv0/2;
    vartmp.h = var.h + dt * th0/2;
    
    [tu1, tv1, th1, vartmp]=tendencies(vartmp, grd, par);
    
    t2 = t0 + dt/2;
    vartmp.u = var.u + dt * tu1/2;
    vartmp.v = var.v + dt * tv1/2;
    vartmp.h = var.h + dt * th1/2;
    
    [tu2, tv2, th2, vartmp]=tendencies(vartmp, grd, par);
    
    
    t3 = t0 + dt;
    vartmp.u = var.u + dt * tu2;
    vartmp.v = var.v + dt * tv2;
    vartmp.h = var.h + dt * th2;
    
    [tu3, tv3, th3, varf]=tendencies(vartmp, grd, par);
    
    if par.tc == 7
        varf.h = var.h + dt*( th0 + 2*th1 + 2*th2 + th3 )/6;
        varf.u = var.u ; %+ dt*( tu0 + 2*tu1 + 2*tu2 + tu3 )/6;
        varf.v = var.v ; %+ dt*( tv0 + 2*tv1 + 2*tv2 + tv3 )/6;
    else
        varf.h = var.h + dt*( th0 + 2*th1 + 2*th2 + th3 )/6;
        varf.u = var.u + dt*( tu0 + 2*tu1 + 2*tu2 + tu3 )/6;
        varf.v = var.v + dt*( tv0 + 2*tv1 + 2*tv2 + tv3 )/6;
        
    end
    
    %  varf.h = var.h + dt*( th0 + 2*th1 + 2*th2 + th3 )/6;
    
end

end