%% Calculate h at u, v and q points
function [hu, hv, hq]=calc_h_uvq(h, grd, mtd)
hu=h;
hv=h;
hq=h;

% bar x of h to u points
hu(2:grd.nx,1:grd.ny)=(h(1:grd.nx-1,1:grd.ny) + h(2:grd.nx,1:grd.ny))/2;
hu(1       ,1:grd.ny)=(h(grd.nx    ,1:grd.ny) + h(1       ,1:grd.ny))/2;

% bar y of h to v points
hv(1:grd.nx,2:grd.ny)=(h(1:grd.nx,1:grd.ny-1) + h(1:grd.nx,2:grd.ny))/2;
hv(1:grd.nx,1       )=(h(1:grd.nx, grd.ny   ) + h(1:grd.nx, 1      ))/2;

% bar xy of h to q points (equivalent to bar y of hu to q points)
hq(1:grd.nx,2:grd.ny)=(hu(1:grd.nx,1:grd.ny-1) + hu(1:grd.nx,2:grd.ny))/2;
hq(1:grd.nx,1       )=(hu(1:grd.nx, grd.ny   ) + hu(1:grd.nx, 1      ))/2;


% Correct scheme if modified method to be used
if mtd == 1
    
   % hu correction ----------------  
   h2=h;

   %bar y of hq back to u points
   h2(1:grd.nx,1:grd.ny-1)=(hq(1:grd.nx,2:grd.ny) + hq(1:grd.nx,1:grd.ny-1))/2;
   h2(1:grd.nx, grd.ny   )=(hq(1:grd.nx, 1      ) + hq(1:grd.nx, grd.ny   ))/2;
   
   %Correct scheme for hu
   hu=hu/3+2*h2/3;
   
   % hv correction ---------------

   %bar x of hq back to v points
   h2(1:grd.nx-1,1:grd.ny)=(hq(2:grd.nx,1:grd.ny) + hq(1:grd.nx-1,1:grd.ny))/2;
   h2(    grd.nx,1:grd.ny)=(hq(1       ,1:grd.ny) + hq( grd.nx   ,1:grd.ny))/2;
   
   %Correct scheme fo hv
   hv=hv/3+2*h2/3;
   
end

% hq(2:grd.nx,2:grd.ny)=(h(1:grd.nx-1,1:grd.ny-1) + h(1:grd.nx-1,2:grd.ny)+h(2:grd.nx,1:grd.ny-1) + h(2:grd.nx,2:grd.ny))/4;
% hq(1       ,2:grd.ny)=(h(grd.nx    ,1:grd.ny-1) + h(grd.nx    ,2:grd.ny)+h(1       ,1:grd.ny-1) + h(1       ,2:grd.ny))/4;
% hq(2:grd.nx,1       )=(h(1:grd.nx-1,    grd.ny) + h(1:grd.nx-1, 1      )+h(2:grd.nx, grd.ny   ) + h(2:grd.nx,  1     ))/4;
% hq(1       ,1       )=(h( grd.nx   ,    grd.ny) + h( grd.nx   , 1      )+h(1       , grd.ny   ) + h(1       ,  1     ))/4;

% %Mid of domain
% for ix=2:grd.nx-1
%     %ixp1=modn(ix+1, grd.nx);
%     ixm1=modn(ix-1, grd.nx);
%     
%     for iy=2:grd.ny-1
%         %iyp1=modn(iy+1, grd.ny);
%         iym1=modn(iy-1, grd.ny);
%         
%         hu(ix,iy)=(h(ixm1, iy) + h(ix, iy))/2;
%         hv(ix,iy)=(h(ix, iym1) + h(ix, iy))/2;
%         hq(ix,iy)=(h(ixm1, iym1) + h(ixm1, iy)+h(ix, iym1) + h(ix, iy))/4;
%     end
% end
% 
% % Borders
% for k=1:grd.nb
%     ix=grd.ixb(k);
%     iy=grd.iyb(k);
%     
%     %ixp1=modn(ix+1, grd.nx);
%     ixm1=modn(ix-1, grd.nx);
%         
%     %iyp1=modn(iy+1, grd.ny);
%     iym1=modn(iy-1, grd.ny);
%     
%     hu(ix,iy)=(h(ixm1, iy) + h(ix, iy))/2;
%     hv(ix,iy)=(h(ix, iym1) + h(ix, iy))/2;
%     hq(ix,iy)=(h(ixm1, iym1) + h(ixm1, iy)+h(ix, iym1) + h(ix, iy))/4;
%     
% end

end