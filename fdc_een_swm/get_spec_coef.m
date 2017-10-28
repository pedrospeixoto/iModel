function w=get_spec_coef(ikappa, jlambda, f, f0, pert, grd)

nxy=grd.nx*grd.ny;

s=fft2((f-f0)/pert)/nxy;
figure(5)
contourf(real(s))
colorbar
figure(6)
contourf(imag(s))
colorbar

if ikappa==0
   w=s(ikappa+1, jlambda+1)+s(ikappa+1, grd.ny-jlambda+1);
elseif jlambda==0
   w=s(ikappa+1, jlambda+1)+s(grd.nx-ikappa+1, jlambda+1);
else
   w=s(ikappa+1, jlambda+1)+s(grd.nx-ikappa+1, grd.ny-jlambda+1);
end

end