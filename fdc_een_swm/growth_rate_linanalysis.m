
function gr=growth_rate_linanalysis(ikappa, jlambda, grd, par, var0, tu0, tv0, th0)

%ikappa=2;
%jlambda=10;


var=var0;
nxy=(grd.ny*grd.nx);
pert=0.000001;

%ikappa=1;
%jlambda=1;
%kappa=2*pi*ikappa/grd.nx;
%lambda=2*pi*jlambda/grd.ny;
mode(1:grd.nx, 1:grd.ny)=0.;
for i=1:grd.nx
    for j=1:grd.ny
        x=(i-1)*grd.dx;
        y=(j-1)*grd.dy;
        kappa=2*pi*ikappa/grd.nx;
        lambda=2*pi*jlambda/grd.ny;
        mode(i, j)=(exp(1i*(kappa*x/grd.dx+lambda*y/grd.dy)));
        %mode(i, j)=mode(i,j) + real(exp(2*pi*1i*((grd.nx-ikappa)*x+(jlambda)*y)));
        %mode(i, j)=mode(i,j) + exp(2*pi*1i*(ikappa*x+(grd.ny-jlambda)*y));
        %mode(i, j)=mode(i,j) + exp(2*pi*1i*((grd.nx-ikappa)*x+(grd.ny-jlambda)*y));
    end
end

% u velocity
var.u=var0.u+pert*mode;
var.v=var0.v;
var.h=var0.h;
[tu, tv, th, var2]=tendencies(var, grd, par);

%s=fft2((tu-tu0)/pert)/nxy;

%w=get_spec_coef(ikappa, jlambda, tu, tu0, pert, grd);
%w=s(ikappa+1, jlambda+1); %+s(grd.nx-ikappa+1, grd.ny-jlambda+1);
wu=fft2((tu-tu0)/pert)/nxy;
wv=fft2((tv-tv0)/pert)/nxy;
wh=fft2((th-th0)/pert)/nxy;

M(1:3,1:3)=0;
M(1,1)=wu(ikappa+1, jlambda+1);
M(2,1)=wv(ikappa+1, jlambda+1);
M(3,1)=wh(ikappa+1, jlambda+1);

% v velocity
var.h=var0.h;
var.u=var0.u;
var.v=var0.v+pert*mode;
[tu, tv, th, var]=tendencies(var, grd, par);

wu=fft2((tu-tu0)/pert)/nxy;
wv=fft2((tv-tv0)/pert)/nxy;
wh=fft2((th-th0)/pert)/nxy;

M(1,2)=wu(ikappa+1, jlambda+1);
M(2,2)=wv(ikappa+1, jlambda+1);
M(3,2)=wh(ikappa+1, jlambda+1);


% h
var.h=var0.h+pert*mode;
var.u=var0.u;
var.v=var0.v;
[tu, tv, th, var]=tendencies(var, grd, par);

wu=fft2((tu-tu0)/pert)/nxy;
wv=fft2((tv-tv0)/pert)/nxy;
wh=fft2((th-th0)/pert)/nxy;

M(1,3)=wu(ikappa+1, jlambda+1);
M(2,3)=wv(ikappa+1, jlambda+1);
M(3,3)=wh(ikappa+1, jlambda+1);

M=M/par.f0;

%Growth rates
gr=max(abs(real(eig(M))));



end


