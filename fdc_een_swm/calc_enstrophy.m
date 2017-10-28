%% Calculate the kinetic energy
function ten=calc_enstrophy(var, grd, par)
h=var.h;
% bar x of h to u points
hu(2:grd.nx,1:grd.ny)=(h(1:grd.nx-1,1:grd.ny) + h(2:grd.nx,1:grd.ny))/2;
hu(1       ,1:grd.ny)=(h(grd.nx    ,1:grd.ny) + h(1       ,1:grd.ny))/2;

% bar xy of h to q points (equivalent to bar y of hu to q points)
hq(1:grd.nx,2:grd.ny)=(hu(1:grd.nx,1:grd.ny-1) + hu(1:grd.nx,2:grd.ny))/2;
hq(1:grd.nx,1       )=(hu(1:grd.nx, grd.ny   ) + hu(1:grd.nx, 1      ))/2;

%Pot enstrophy
ten=sum(sum(var.q.*var.q.*hq))*grd.dx*grd.dy;
% Total enstrophy
%ten=sum(sum(var.q.*var.q))*grd.dx*grd.dy;

end