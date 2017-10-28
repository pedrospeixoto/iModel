%% Calculate the kinetic energy - unused
function gradmax=calc_grad_kenergy(var, grd, par)

gradx(1:grd.nx-1,1:grd.ny)=(var.ke(2:grd.nx,1:grd.ny)-var.ke(1:grd.nx-1,1:grd.ny))/grd.dx;
gradx( grd.nx   ,1:grd.ny)=(var.ke(  1     ,1:grd.ny)-var.ke( grd.nx   ,1:grd.ny))/grd.dx;

divy(1:grd.nx,1:grd.ny-1)=(var.ke(1:grd.nx,2:grd.ny)-var.ke(1:grd.nx,1:grd.ny-1))/grd.dy;
divy(1:grd.nx, grd.ny   )=(var.ke(1:grd.nx,1       )-var.ke(1:grd.nx, grd.ny   ))/grd.dy;

gradke=var.ke

gradmax=max(max(gradke));

end