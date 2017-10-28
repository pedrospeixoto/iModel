
%% Contourf variables
% Select variable to be ploted using string in 'field'
function plot_vars_t(t, var, grd, par, fig)

%fig=figure('Color',[1 1 1]);
figure(fig)
cla;
subplot(2,2,1)
contourf(grd.xh, grd.yh, var.h', 20, 'LineColor','none')
ylabel('y','fontsize',13)
xlabel('x','fontsize',13)
axis([0,1,0,1])
title([ 'h  t = ', num2str(t,'%8.4f')], 'fontsize',12)
colorbar

subplot(2,2,2)        
contourf(grd.xz, grd.yz, var.z', 20, 'LineColor','none')
ylabel('y','fontsize',13)
xlabel('x','fontsize',13)
axis([0,1,0,1])
title([ 'z  t = ', num2str(t,'%8.4f')], 'fontsize',12)
colorbar

subplot(2,2,3)        
contourf(grd.xu, grd.yu, var.u', 20, 'LineColor','none')
ylabel('y','fontsize',13)
xlabel('x','fontsize',13)
axis([0,1,0,1])
title([ 'u  t = ', num2str(t,'%8.4f')], 'fontsize',12)
colorbar

subplot(2,2,4)        
    
contourf(grd.xv, grd.yv, var.v', 20, 'LineColor','none')
%Graph labels
ylabel('y','fontsize',13)
xlabel('x','fontsize',13)
axis([0,1,0,1])
title(['v  t = ', num2str(t,'%8.4f')], 'fontsize',12)
colorbar

drawnow
% coment=annotation('textbox',[0.55 0.01 1 0.05]);
% set(coment,'linestyle','none')
% set(coment,'FitBoxToText','on')
% set(coment,'string',[...
%     '  t = ', num2str(t), ...
%     '  f_0 = ', num2str(par.f0), ...
%     '  n_x = ', num2str(grd.nx) ...
%     ]);

end