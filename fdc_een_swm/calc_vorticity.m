%% Calculate Vorticity
function z=calc_vorticity(u, v, grd)

z(1:grd.nx,1:grd.ny)=0;
z(2:grd.nx,2:grd.ny)=(v(2:grd.nx,2:grd.ny) - v(1:grd.nx-1,2:grd.ny))/grd.dx - (u(2:grd.nx,2:grd.ny) - u(2:grd.nx,1:grd.ny-1))/grd.dy;

%Borders
z(1,2:grd.ny)=(v(1,2:grd.ny) - v(grd.nx,2:grd.ny))/grd.dx - (u(1,2:grd.ny) - u(1,1:grd.ny-1))/grd.dy;
z(2:grd.nx,1)=(v(2:grd.nx,1) - v(1:grd.nx-1,1))/grd.dx - (u(2:grd.nx,1) - u(2:grd.nx,grd.ny))/grd.dy;
z(1,1)=(v(1,1) - v(grd.nx,1))/grd.dx - (u(1,1) - u(1,grd.ny))/grd.dy;

% %Mid of domain
% for ix=2:grd.nx-1
%     ixm1=ix-1; 
%     %ixm1=m1mod(ix,grd.nx);
%     %ixm1=modn(ix-1, grd.nx);
%     ixp1=ix;
%     for iy=2:grd.ny-1
%         %iym1= m1mod(iy,grd.ny); 
%         %iym1= modn(iy-1, grd.ny);
%         iym1= iy-1;
%         iyp1=iy;
%         
%         z(ix,iy)=(v(ixp1, iy) - v(ixm1, iy))/grd.dx - (u(ix, iyp1) - u(ix, iym1))/grd.dy;
%     end
% end

% %Border
% for k=1:grd.nb
%     ix=grd.ixb(k);
%     iy=grd.iyb(k);
%     
%     ixm1=modn(ix-1, grd.nx);
%     ixp1=ix;
%     iym1= modn(iy-1, grd.ny);
%     iyp1=iy;
%     
%     z(ix,iy)=(v(ixp1, iy) - v(ixm1, iy))/grd.dx - (u(ix, iyp1) - u(ix, iym1))/grd.dy;
% end

% Full domain - slow
% for ix=1:grd.nx
%     ixm1=m1mod(ix,grd.nx);
%     %ixm1=modn(ix-1, grd.nx);
%     ixp1=ix;
%     for iy=1:grd.ny
%         iym1= m1mod(iy,grd.ny); 
%         %iym1= modn(iy-1, grd.ny);
%         iyp1=iy;
%         
%         z(ix,iy)=(v(ixp1, iy) - v(ixm1, iy))/grd.dx - (u(ix, iyp1) - u(ix, iym1))/grd.dy;
%     end
% end

% Separate all border parts - ugly
% % (0,0) corner
% ix=1;
% ixm1=grd.nx;
% iy=1;
% iym1=grd.ny;
% iyp1=iy;
% z(ix,iy)=(v(ixp1, iy) - v(ixm1, iy))/grd.dx - (u(ix, iyp1) - u(ix, iym1))/grd.dy;
% 
% % (0,y) border
% 
% ix=1;
% ixm1=grd.nx;
% %ixm1=m1mod(ix,grd.nx);
% %ixm1=modn(ix-1, grd.nx);
% ixp1=ix;
% for iy=2:grd.ny
%     %iym1= m1mod(iy,grd.ny);
%     %iym1= modn(iy-1, grd.ny);
%     iym1= iy-1;
%     iyp1=iy;
%     
%     z(ix,iy)=(v(ixp1, iy) - v(ixm1, iy))/grd.dx - (u(ix, iyp1) - u(ix, iym1))/grd.dy;
% end
% 
% % (x,0) border
% iy=1;
% iym1=grd.ny;
% iyp1=iy;
% for ix=2:grd.nx
%     ixm1= ix-1;
%     ixp1=ix;
%     z(ix,iy)=(v(ixp1, iy) - v(ixm1, iy))/grd.dx - (u(ix, iyp1) - u(ix, iym1))/grd.dy;
% end

end