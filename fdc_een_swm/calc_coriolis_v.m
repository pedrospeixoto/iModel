%% Calculate the Coriolis term tendencies
function zu=calc_coriolis_v(q, qy, qx, us, grd, mtd)

if mtd<2 %Use een scheme
% First term

% at u points
qyus=qy.*us;

% qyus bar y - at q points
qyusy(1:grd.nx,2:grd.ny)=(qyus(1:grd.nx,1:grd.ny-1)+qyus(1:grd.nx,2:grd.ny))/2;
qyusy(1:grd.nx, 1      )=(qyus(1:grd.nx, grd.ny   )+qyus(1:grd.nx, 1      ))/2;

% qyusy bar x - at v points
qyusyx(1:grd.nx-1,1:grd.ny)=(qyusy(2:grd.nx,1:grd.ny)+qyusy(1:grd.nx-1,1:grd.ny))/2;
qyusyx( grd.nx   ,1:grd.ny)=(qyusy(1       ,1:grd.ny)+qyusy( grd.nx   ,1:grd.ny))/2;

%Second term

% us bar y = usy - at q points
usy(1:grd.nx,2:grd.ny)=(us(1:grd.nx,1:grd.ny-1)+us(1:grd.nx,2:grd.ny))/2;
usy(1:grd.nx, 1      )=(us(1:grd.nx, grd.ny   )+us(1:grd.nx, 1      ))/2;

% usy bar x - at v points
usyx(1:grd.nx-1,1:grd.ny)=(usy(2:grd.nx,1:grd.ny)+usy(1:grd.nx-1,1:grd.ny))/2;
usyx( grd.nx   ,1:grd.ny)=(usy(1       ,1:grd.ny)+usy( grd.nx   ,1:grd.ny))/2;


% Third term
qusy=q.*usy;

% Bar x in qusy - at v points
qusyx(1:grd.nx-1,1:grd.ny)=(qusy(2:grd.nx,1:grd.ny)+qusy(1:grd.nx-1,1:grd.ny))/2;
qusyx( grd.nx   ,1:grd.ny)=(qusy(1       ,1:grd.ny)+qusy( grd.nx   ,1:grd.ny))/2;

%Result
zu=(2/3)*qyusyx+(2/3)*qx.*usyx-(1/3)*qusyx;

elseif mtd ==2 % use Energy conserving sadourny scheme

% us bar y = usy - at q points
usy(1:grd.nx,2:grd.ny)=(us(1:grd.nx,1:grd.ny-1)+us(1:grd.nx,2:grd.ny))/2;
usy(1:grd.nx, 1      )=(us(1:grd.nx, grd.ny   )+us(1:grd.nx, 1      ))/2;
% Third term
qusy=q.*usy;

% Bar x in qusy - at v points
qusyx(1:grd.nx-1,1:grd.ny)=(qusy(2:grd.nx,1:grd.ny)+qusy(1:grd.nx-1,1:grd.ny))/2;
qusyx( grd.nx   ,1:grd.ny)=(qusy(1       ,1:grd.ny)+qusy( grd.nx   ,1:grd.ny))/2;

%Result
zu=qusyx;

elseif mtd==3 %Enstrophy conserving sadourny scheme
    
% us bar y = usy - at q points
usy(1:grd.nx,2:grd.ny)=(us(1:grd.nx,1:grd.ny-1)+us(1:grd.nx,2:grd.ny))/2;
usy(1:grd.nx, 1      )=(us(1:grd.nx, grd.ny   )+us(1:grd.nx, 1      ))/2;

% usy bar x - at v points
usyx(1:grd.nx-1,1:grd.ny)=(usy(2:grd.nx,1:grd.ny)+usy(1:grd.nx-1,1:grd.ny))/2;
usyx( grd.nx   ,1:grd.ny)=(usy(1       ,1:grd.ny)+usy( grd.nx   ,1:grd.ny))/2;

zu=qx.*usyx;


end
% %Initialize variables
% qyus=us;
% qxusbarxy=us;
% qusy=us;
% qyusbarxy=us;
% qusybary=us;
% 
% for ix=1:grd.nx
%     ixp1=modn(ix+1, grd.nx);
%     %ixm1=modn(ix-1, grd.nx);
%     
%     for iy=1:grd.ny
%         
%         iyp1=modn(iy+1, grd.ny);
%         iym1=modn(iy-1, grd.ny);
%         qyus0(ix,iy)=us(ix,iy)*(q(ix, iyp1)+q(ix,iy))/2;
%         qxusbarxy0(ix,iy)=(2/3)*((q(ix, iy)+q(ixp1,iy))/2)*(us(ix, iy)+us(ix,iym1)+us(ixp1, iym1)+us(ixp1,iy))/4;
%         qusy0(ix,iy)=q(ix,iy)*(us(ix, iy)+us(ix, iym1))/2;
%     end
% end
% 
% qxusbarxy0-(2/3)*qx.*usyx

% 
% % Do boundaries
% for k=1:grd.nb
%     ix=grd.ixb(k);
%     iy=grd.iyb(k);
%     
%     ixp1=modn(ix+1, grd.nx);
%     %ixm1=modn(ix-1, grd.nx);
%     
%     iyp1=modn(iy+1, grd.ny);
%     iym1=modn(iy-1, grd.ny);
%     
%     qyus(ix,iy)=us(ix,iy)*(q(ix, iyp1)+q(ix,iy))/2;
%     qxusbarxy(ix,iy)=(2/3)*((q(ix, iy)+q(ixp1,iy))/2)*(us(ix, iy)+us(ix,iym1)+us(ixp1, iym1)+us(ixp1,iy))/4;
%     qusy(ix,iy)=q(ix,iy)*(us(ix, iy)+us(ix, iym1))/2;
%     
% end

% % Mid domain
% for ix=2:grd.nx-1
%     ixp1=modn(ix+1, grd.nx);
%     %ixm1=modn(ix-1, grd.nx);
%     
%     for iy=2:grd.ny-1
%         %iyp1=modn(iy+1, grd.ny);
%         iym1=modn(iy-1, grd.ny);
%         qyusbarxy(ix,iy)=(2/3)*(qyus(ix, iy)+qyus(ix,iym1)+qyus(ixp1, iym1)+qyus(ixp1,iy))/4;
%         qusybary(ix,iy)=(-1/3)*(qusy(ix, iy)+qusy(ixp1, iy))/2;
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
%     %iyp1=modn(iy+1, grd.ny);
%     iym1=modn(iy-1, grd.ny);
%     
%     qyusbarxy(ix,iy)=(2/3)*(qyus(ix, iy)+qyus(ix,iym1)+qyus(ixp1, iym1)+qyus(ixp1,iy))/4;
%     qusybary(ix,iy)=(-1/3)*(qusy(ix, iy)+qusy(ixp1, iy))/2;
%     
% end
% 
% %Result
% zu=qyusbarxy+qxusbarxy+qusybary;
end