%% Calculate the gradient of the Bernoulli potential (g(h+b)+K)
function gy=calc_grady(phi, grd)
%grav    = 9.80616 ; 10;

gy=phi;
gy(1:grd.nx,2:grd.ny)=phi(1:grd.nx,2:grd.ny)- phi(1:grd.nx,1:grd.ny-1);
gy(1:grd.nx,1)=phi(1:grd.nx,1)- phi(1:grd.nx,grd.ny);

gy=gy./grd.dy;

% % Mid domain
% for ix=2:grd.nx-1
%     %ixp1=modn(ix+1, grd.nx);
%     %ixm1=modn(ix-1, grd.nx);
%     
%     for iy=2:grd.ny-1
%         %iyp1=modn(iy+1, grd.ny);
%         iym1=iy-1; %modn(iy-1, grd.ny);
%         gy(ix,iy)=(phi(ix,iy)- phi(ix,iym1));
%     end
% end
% 
% %Border
% for k=1:grd.nb
%     ix=grd.ixb(k);
%     iy=grd.iyb(k);
%     
%     %ixp1=modn(ix+1, grd.nx);
%     %ixm1=modn(ix-1, grd.nx);
%         
%     %iyp1=modn(iy+1, grd.ny);
%     iym1=modn(iy-1, grd.ny);
%     
%     gy(ix,iy)=(phi(ix,iy)- phi(ix,iym1));
%         
% end



end