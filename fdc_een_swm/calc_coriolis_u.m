%% Calculate the Coriolis term tendencies
function zv=calc_coriolis_u(q, qy, qx, vs, grd, mtd)

if mtd<2 %Use een scheme
% First term

% at v points
qxvs=qx.*vs;

% qxvs bar x - at q points
qxvsx(2:grd.nx,1:grd.ny)=(qxvs(2:grd.nx,1:grd.ny)+qxvs(1:grd.nx-1,1:grd.ny))/2;
qxvsx( 1      ,1:grd.ny)=(qxvs(1       ,1:grd.ny)+qxvs( grd.nx   ,1:grd.ny))/2;

% qxvsx bar y - at u points
qxvsxy(1:grd.nx,1:grd.ny-1)=(qxvsx(1:grd.nx,1:grd.ny-1)+qxvsx(1:grd.nx,2:grd.ny))/2;
qxvsxy(1:grd.nx, grd.ny   )=(qxvsx(1:grd.nx, grd.ny   )+qxvsx(1:grd.nx, 1      ))/2;

% Second term

% vs bar x - at q points
vsx(2:grd.nx,1:grd.ny)=(vs(2:grd.nx,1:grd.ny)+vs(1:grd.nx-1,1:grd.ny))/2;
vsx( 1      ,1:grd.ny)=(vs(1       ,1:grd.ny)+vs( grd.nx   ,1:grd.ny))/2;

% vsx bar y = vsxy - at u points
vsxy(1:grd.nx,1:grd.ny-1)=(vsx(1:grd.nx,1:grd.ny-1)+vsx(1:grd.nx,2:grd.ny))/2;
vsxy(1:grd.nx, grd.ny   )=(vsx(1:grd.nx, grd.ny   )+vsx(1:grd.nx, 1      ))/2;

% Third term
qvsx=q.*vsx;

qvsxy(1:grd.nx,1:grd.ny-1)=(qvsx(1:grd.nx,1:grd.ny-1)+qvsx(1:grd.nx,2:grd.ny))/2;
qvsxy(1:grd.nx, grd.ny   )=(qvsx(1:grd.nx, grd.ny   )+qvsx(1:grd.nx, 1      ))/2;

%Result
zv=(2/3)*qxvsxy+(2/3)*qy.*vsxy-(1/3)*qvsxy;

elseif mtd ==2 % use Energy conserving sadourny scheme

% vs bar x - at q points
vsx(2:grd.nx,1:grd.ny)=(vs(2:grd.nx,1:grd.ny)+vs(1:grd.nx-1,1:grd.ny))/2;
vsx( 1      ,1:grd.ny)=(vs(1       ,1:grd.ny)+vs( grd.nx   ,1:grd.ny))/2;

% Third term
qvsx=q.*vsx;

qvsxy(1:grd.nx,1:grd.ny-1)=(qvsx(1:grd.nx,1:grd.ny-1)+qvsx(1:grd.nx,2:grd.ny))/2;
qvsxy(1:grd.nx, grd.ny   )=(qvsx(1:grd.nx, grd.ny   )+qvsx(1:grd.nx, 1      ))/2;

%Result
zv=qvsxy;

elseif mtd==3 %Enstrophy conserving sadourny scheme
    
    
% vs bar x - at q points
vsx(2:grd.nx,1:grd.ny)=(vs(2:grd.nx,1:grd.ny)+vs(1:grd.nx-1,1:grd.ny))/2;
vsx( 1      ,1:grd.ny)=(vs(1       ,1:grd.ny)+vs( grd.nx   ,1:grd.ny))/2;

% vsx bar y = vsxy - at u points
vsxy(1:grd.nx,1:grd.ny-1)=(vsx(1:grd.nx,1:grd.ny-1)+vsx(1:grd.nx,2:grd.ny))/2;
vsxy(1:grd.nx, grd.ny   )=(vsx(1:grd.nx, grd.ny   )+vsx(1:grd.nx, 1      ))/2;

zv=qy.*vsxy;

end
% % Do the mid-domain calculation
% 
% %Calculate qbarx * vs
% %Calculate qbary * vsbarxy
% %Calculate q * vsbarx
% for ix=2:grd.nx-1
%     ixp1=ix+1; %modn(ix+1, grd.nx);
%     ixm1=ix-1; %modn(ix-1, grd.nx);
%     
%     for iy=2:grd.ny-1
%         iyp1=iy+1; %modn(iy+1, grd.ny);
%         %iym1=modn(iy-1, grd.ny);
%         qxvs(ix,iy)=vs(ix,iy)*(q(ixp1, iy)+q(ix,iy))/2;
%         qyvsbarxy(ix,iy)=(2/3)*((q(ix, iyp1)+q(ix,iy))/2)*(vs(ixm1, iy)+vs(ixm1,iyp1)+vs(ix, iy)+vs(ix,iyp1))/4;
%         qvsx(ix,iy)=q(ix,iy)*(vs(ixm1, iy)+vs(ix, iy))/2;
%     end
% end
% 
% % Calculate border terms
% 
% %Calculate qbarx * vs
% %Calculate qbary * vsbarxy
% %Calculate q * vsbarx
% for k=1:grd.nb
%     ix=grd.ixb(k);
%     iy=grd.iyb(k);
%     
%     ixp1=modn(ix+1, grd.nx);
%     ixm1=modn(ix-1, grd.nx);
%         
%     iyp1=modn(iy+1, grd.ny);
%     %iym1=modn(iy-1, grd.ny);
%     
%     qxvs(ix,iy)=vs(ix,iy)*(q(ixp1, iy)+q(ix,iy))/2;
%     qyvsbarxy(ix,iy)=(2/3)*((q(ix, iyp1)+q(ix,iy))/2)*(vs(ixm1, iy)+vs(ixm1,iyp1)+vs(ix, iy)+vs(ix,iyp1))/4;
%     qvsx(ix,iy)=q(ix,iy)*(vs(ixm1, iy)+vs(ix, iy))/2;
% end
% 
% % Mid domain calculation
% 
% %Calculate (qbarx * vs )barxy
% %Calculate (q * vsx)bary
% for ix=2:grd.nx-1
%     %ixp1=modn(ix+1, grd.nx);
%     ixm1=ix-1; %modn(ix-1, grd.nx);
%     
%     for iy=2:grd.ny-1
%         iyp1=iy+1; %modn(iy+1, grd.ny);
%         %iym1=modn(iy-1, grd.ny);
%         qxvsbarxy(ix,iy)=(2/3)*(qxvs(ixm1, iy)+qxvs(ixm1,iyp1)+qxvs(ix, iy)+qxvs(ix,iyp1))/4;
%         qvsxbary(ix,iy)=(-1/3)*(qvsx(ix, iy)+qvsx(ix, iyp1))/2;
%     end
% end
% 
% % Border
% 
% %Calculate (qbarx * vs )barxy
% %Calculate (q * vsx)bary
% for k=1:grd.nb
%     ix=grd.ixb(k);
%     iy=grd.iyb(k);
%     
%     %ixp1=modn(ix+1, grd.nx);
%     ixm1=modn(ix-1, grd.nx);
%         
%     iyp1=modn(iy+1, grd.ny);
%     %iym1=modn(iy-1, grd.ny);
%     
%     qxvsbarxy(ix,iy)=(2/3)*(qxvs(ixm1, iy)+qxvs(ixm1,iyp1)+qxvs(ix, iy)+qxvs(ix,iyp1))/4;
%     qvsxbary(ix,iy)=(-1/3)*(qvsx(ix, iy)+qvsx(ix, iyp1))/2;
%     
% end
% 
% % Result
% 
% zv=qxvsbarxy+qyvsbarxy+qvsxbary;
end