%% Calculate the divergence term
function div=calc_div(us, vs, grd) 

divx(1:grd.nx-1,1:grd.ny)=(us(2:grd.nx,1:grd.ny)-us(1:grd.nx-1,1:grd.ny))/grd.dx;
divx( grd.nx   ,1:grd.ny)=(us(  1     ,1:grd.ny)-us( grd.nx   ,1:grd.ny))/grd.dx;

divy(1:grd.nx,1:grd.ny-1)=(vs(1:grd.nx,2:grd.ny)-vs(1:grd.nx,1:grd.ny-1))/grd.dy;
divy(1:grd.nx, grd.ny   )=(vs(1:grd.nx,1       )-vs(1:grd.nx, grd.ny   ))/grd.dy;

div=divx+divy;

% 
% %Mid domain
% for ix=2:grd.nx-1
%     ixp1=modn(ix+1, grd.nx);
%     %ixm1=modn(ix-1, grd.nx);
%     
%     for iy=2:grd.ny-1
%         iyp1=modn(iy+1, grd.ny);
%         %iym1=modn(iy-1, grd.ny);
%         
%         div(ix,iy)=(us(ixp1,iy)-us(ix,iy))/grd.dx+(vs(ix,iyp1)-vs(ix,iy))/grd.dy;
%     end
% end
% 
% %Borders
% for k=1:grd.nb
%     ix=grd.ixb(k);
%     iy=grd.iyb(k);
%     
%     ixp1=modn(ix+1, grd.nx);
%     %ixm1=modn(ix-1, grd.nx);
%     
%     iyp1=modn(iy+1, grd.ny);
%     %iym1=modn(iy-1, grd.ny);
%     
%     div(ix,iy)=(us(ixp1,iy)-us(ix,iy))/grd.dx+(vs(ix,iyp1)-vs(ix,iy))/grd.dy;
%     
% end

end