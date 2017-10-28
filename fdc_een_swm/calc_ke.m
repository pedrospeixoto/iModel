%% Calculate the kinetic energy
function ke=calc_ke(u, v, grd, mtd)

u2=u.*u;
v2=v.*v;

u2x(1:grd.nx-1,1:grd.ny)=(u2(1:grd.nx-1,1:grd.ny)+u2(2:grd.nx,1:grd.ny))/2;
u2x( grd.nx   ,1:grd.ny)=(u2(grd.nx    ,1:grd.ny)+u2(   1    ,1:grd.ny))/2;

v2y(1:grd.nx,1:grd.ny-1)=(v2(1:grd.nx,1:grd.ny-1)+v2(1:grd.nx,2:grd.ny))/2;
v2y(1:grd.nx, grd.ny   )=(v2(1:grd.nx, grd.ny   )+v2(1:grd.nx,1       ))/2;

ke=(1/2)*(u2x+v2y);

if mtd == 1 %use corrected scheme (Hollingsworth stable)
    
    % bar y on u2x
    u2xy(1:grd.nx,2:grd.ny)=(u2x(1:grd.nx,2:grd.ny)+u2x(1:grd.nx,1:grd.ny-1))/2;
    u2xy(1:grd.nx, 1      )=(u2x(1:grd.nx, grd.ny )+u2x(1:grd.nx,1       ))/2;
    
    % bar y on u2xy
    u2xyy(1:grd.nx,1:grd.ny-1)=(u2xy(1:grd.nx,2:grd.ny)+u2xy(1:grd.nx,1:grd.ny-1))/2;
    u2xyy(1:grd.nx, grd.ny   )=(u2xy(1:grd.nx, grd.ny )+u2xy(1:grd.nx,1       ))/2;
    
    % bar x on v2y 
    v2yx(2:grd.nx,1:grd.ny)=(v2y(2:grd.nx,1:grd.ny)+v2y(1:grd.nx-1,1:grd.ny))/2;
    v2yx( 1      ,1:grd.ny)=(v2y(grd.nx  ,1:grd.ny)+v2y(   1      ,1:grd.ny))/2;
        
    % bar x on v2yx 
    v2yxx(1:grd.nx-1,1:grd.ny)=(v2yx(2:grd.nx,1:grd.ny)+v2yx(1:grd.nx-1,1:grd.ny))/2;
    v2yxx( grd.nx   ,1:grd.ny)=(v2yx(grd.nx  ,1:grd.ny)+v2yx(   1      ,1:grd.ny))/2;
    
    ke_m=(u2xyy+v2yxx)/2;
    
    %Correct the original scheme
    ke=(1/3)*ke+(2/3)*ke_m;
end

% if mtd == 1 %use corrected scheme (Hollingsworth stable)
%     
%     % First average the u and v to the vorticity points
%     uy2(1:grd.nx,2:grd.ny)=(u(1:grd.nx,2:grd.ny)+u(1:grd.nx,1:grd.ny-1))/2;
%     uy2(1:grd.nx, 1      )=(u(1:grd.nx, grd.ny )+u(1:grd.nx,1       ))/2;
%     uy2=uy2.*uy2;
%         
%     vx2(2:grd.nx,1:grd.ny)=(v(2:grd.nx,1:grd.ny)+v(1:grd.nx-1,1:grd.ny))/2;
%     vx2( 1      ,1:grd.ny)=(v(grd.nx    ,1:grd.ny)+v(   1    ,1:grd.ny))/2;
%     vx2=vx2.*vx2;
%     
%     ke_z=(uy2+vx2)/2;
%     
%     
%     %Now average ke_z to the h point
%     % Start putting into v points (bar x)
%     ke_v(1:grd.nx-1,1:grd.ny)=(ke_z(1:grd.nx-1,1:grd.ny)+ke_z(2:grd.nx,1:grd.ny))/2;
%     ke_v( grd.nx   ,1:grd.ny)=(ke_z(grd.nx    ,1:grd.ny)+ke_z(   1    ,1:grd.ny))/2;
%     
%     %Now average to h points
%     ke_h(1:grd.nx,1:grd.ny-1)=(ke_v(1:grd.nx,1:grd.ny-1)+ke_v(1:grd.nx,2:grd.ny))/2;
%     ke_h(1:grd.nx, grd.ny   )=(ke_v(1:grd.nx, grd.ny   )+ke_v(1:grd.nx,1       ))/2;
%     
%     %Correct the original scheme
%     ke=(1/3)*ke+(2/3)*ke_h;
% end

% %+(v2(ix,iy)+v2(ix,iyp1))/2;
% 
% %Mid domain
% for ix=2:grd.nx-1
%     ixp1=ix+1; %modn(ix+1, grd.nx);
%     %ixm1=modn(ix-1, grd.nx);
%     
%     for iy=2:grd.ny-1
%         iyp1=iy+1; %modn(iy+1, grd.ny);
%         %iym1=modn(iy-1, grd.ny);
%         ke(ix,iy)=(u2(ix,iy)+u2(ixp1,iy))/2+(v2(ix,iy)+v2(ix,iyp1))/2;
%     end
% end
% 
% %Boundaries
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
%     ke(ix,iy)=(u2(ix,iy)+u2(ixp1,iy))/2+(v2(ix,iy)+v2(ix,iyp1))/2;
% end
% 
%ke=ke/2;
% 
% if mtd == 1 %use corrected scheme (Hollingsworth stable)
%     
%     % First average the ke to the vorticity points
%     ke_z=ke;
%     %Mid domain
%     for ix=2:grd.nx-1
%         %ixp1=ix+1; %modn(ix+1, grd.nx);
%         ixm1=ix-1; %modn(ix-1, grd.nx);
%         
%         for iy=2:grd.ny-1
%             %iyp1=iy+1; %modn(iy+1, grd.ny);
%             iym1=iy-1; %modn(iy-1, grd.ny);
%             ke_z(ix,iy)=(1/2)*(((u(ix,iy)+u(ix,iym1))/2)^2+((v(ix,iy)+v(ixm1,iy))/2)^2);
%         end
%     end
%     
%     %Boundaries
%     for k=1:grd.nb
%         ix=grd.ixb(k);
%         iy=grd.iyb(k);
%         
%         %ixp1=modn(ix+1, grd.nx);
%         ixm1=modn(ix-1, grd.nx);
%         
%         %iyp1=modn(iy+1, grd.ny);
%         iym1=modn(iy-1, grd.ny);
%         
%         ke_z(ix,iy)=(1/2)*(((u(ix,iy)+u(ix,iym1))/2)^2+((v(ix,iy)+v(ixm1,iy))/2)^2);
%     end
%     
%     %Now average them to the h point
%     
%     ke_h=ke;
%     %Mid domain
%     for ix=2:grd.nx-1
%         ixp1=ix+1; %modn(ix+1, grd.nx);
%         %ixm1=ix-1; %modn(ix-1, grd.nx);
%         
%         for iy=2:grd.ny-1
%             iyp1=iy+1; %modn(iy+1, grd.ny);
%             %iym1=iy-1; %modn(iy-1, grd.ny);
%             ke_h(ix,iy)=(1/4)*(ke_z(ix,iy)+ke_z(ixp1,iy)+ke_z(ixp1,iyp1)+ke_z(ix,iyp1) );
%         end
%     end
%     
%     %Boundaries
%     for k=1:grd.nb
%         ix=grd.ixb(k);
%         iy=grd.iyb(k);
%         
%         ixp1=modn(ix+1, grd.nx);
%         %ixm1=modn(ix-1, grd.nx);
%         
%         iyp1=modn(iy+1, grd.ny);
%         %iym1=modn(iy-1, grd.ny);
%         
%         ke_h(ix,iy)=(1/4)*(ke_z(ix,iy)+ke_z(ixp1,iy)+ke_z(ixp1,iyp1)+ke_z(ix,iyp1) );
%     end
%     
%     %Correct the original scheme
%     ke=(1/3)*ke+(2/3)*ke_h;
% end



end