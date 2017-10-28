%% Calculate the kinetic energy
function tenergy=calc_energy(var, grd, par)

kin_en=var.ke*grd.dx*grd.dy;; %var.h.*var.ke*grd.dx*grd.dy;
pot_en=par.g*var.h.*((1/2)*var.h+var.b)*grd.dx*grd.dy;
%tenergy=sum(sum(kin_en))+sum(sum(pot_en));
tenergy=sum(sum(kin_en)); %+sum(sum(pot_en));
end