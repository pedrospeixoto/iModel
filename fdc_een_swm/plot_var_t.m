
%% Contourf variables
% Select variable to be ploted using string in 'field'
function plot_var_t(t, varplot, grd, par, pos, plot_title, fig)

%fig=figure('Color',[1 1 1]);
figure(fig)
cla;
switch pos
    case {'h'}     
        contourf(grd.xh, grd.yh, varplot', 100, 'LineColor','none')
    case {'z', 'q'}     
        contourf(grd.xz, grd.yz, varplot', 100, 'LineColor','none')
    case {'u'}     
        contourf(grd.xu, grd.yu, varplot', 100, 'LineColor','none')
    case {'v'}     
        contourf(grd.xv, grd.yv, varplot', 100, 'LineColor','none')
end
%Graph labels
ylabel('y','fontsize',13)
xlabel('x','fontsize',13)
axis([0,1,0,1])
title([plot_title, ' t = ', num2str(t,'%8.4f')], 'fontsize',12)
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