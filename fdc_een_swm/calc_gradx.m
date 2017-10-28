%% Calculate the gradient of the Bernoulli potential (g(h+b)+K)
function gx=calc_gradx(phi, grd)

gx=phi;
gx(2:grd.nx,1:grd.ny)=phi(2:grd.nx,1:grd.ny)- phi(1:grd.nx-1,1:grd.ny);
gx(1,1:grd.ny)=phi(1,1:grd.ny)- phi(grd.nx,1:grd.ny);
gx=gx./grd.dx;

% % Mid domain
% for ix=2:grd.nx-1
%     %ixp1=modn(ix+1, grd.nx);
%     ixm1=ix-1; %modn(ix-1, grd.nx);
%     
%     for iy=2:grd.ny-1
%         %iyp1=modn(iy+1, grd.ny);
%         %iym1=modn(iy-1, grd.ny);
%         gx(ix,iy)=(phi(ix,iy)- phi(ixm1,iy));
%     end
% end
% 
% 
% %Borders
% for k=1:grd.nb
%     ix=grd.ixb(k);
%     iy=grd.iyb(k);
%     
%     %ixp1=modn(ix+1, grd.nx);
%     ixm1=modn(ix-1, grd.nx);
%     
%     %iyp1=modn(iy+1, grd.ny);
%     %iym1=modn(iy-1, grd.ny);
%     gx(ix,iy)=(phi(ix,iy)- phi(ixm1,iy));
%     
% end
% 
% 


end