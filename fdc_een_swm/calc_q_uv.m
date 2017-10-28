%% Calculate h at u, v and q points
function [qu, qv]=calc_q_uv(q, grd)
qu=q;
qv=q;

% q bar y - at u points
qu(1:grd.nx,1:grd.ny-1)=(q(1:grd.nx,1:grd.ny-1)+q(1:grd.nx,2:grd.ny))/2;
qu(1:grd.nx, grd.ny   )=(q(1:grd.nx, grd.ny   )+q(1:grd.nx, 1      ))/2;


% q bar x - at v points
qv(1:grd.nx-1,1:grd.ny)=(q(2:grd.nx,1:grd.ny)+q(1:grd.nx-1,1:grd.ny))/2;
qv(  grd.nx  ,1:grd.ny)=(q(1       ,1:grd.ny)+q( grd.nx   ,1:grd.ny))/2;

end