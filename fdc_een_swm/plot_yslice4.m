%% Plot y slice
function varmax=plot_yslice4(t, k, var, var0, varmax, grd, fig)

        figure(fig)

        subplot(2,2,1)
        error=var.h(grd.nx/2,:)-var0.h(grd.nx/2,:);
        for i=1:grd.ny
            if (error(i)>varmax.hp(i)) && varmax.hp(i) >= 0
                varmax.hp(i)=error(i);
            elseif (error(i)<varmax.hn(i)) && varmax.hn(i) <= 0
                 varmax.hn(i)=error(i);
            end
        end
        %plot(grd.yh, error+var0.h(grd.nx/2,:))
        %Sum the initial h
        %plot(grd.yh, error+var0.h(grd.nx/2,:), grd.yh, ...
        %    varmax.hn+var0.h(grd.nx/2,:), grd.yh, varmax.hp+var0.h(grd.nx/2,:))
        %Plot just error
        plot(grd.yh, error, grd.yh, ...
            varmax.hn, grd.yh, varmax.hp)
        title('h ')
        %xlabel('y (x=0.5)')
        
        subplot(2,2,2)
        error=var.u(grd.nx/2,:)-var0.u(grd.nx/2,:);
        for i=1:grd.ny
            if (error(i)>varmax.up(i)) && varmax.up(i) >= 0
                varmax.up(i)=error(i);
            elseif (error(i)<varmax.un(i)) && varmax.un(i) <= 0
                 varmax.un(i)=error(i);
            end
        end
        %plot(grd.yu, error)
        plot(grd.yu, error, grd.yu, varmax.un, grd.yu, varmax.up)
        plot(grd.yu, var.u(grd.nx/2,:)-var0.u(grd.nx/2,:))
        title('u')
        %xlabel('y (x=0.5)')
        
        subplot(2,2,3)
        for i=1:grd.ny
            if (var.v(grd.nx/2,i)>varmax.vp(i)) && varmax.vp(i) >= 0
                varmax.vp(i)=var.v(grd.nx/2,i);
            elseif (var.v(grd.nx/2,i)<varmax.vn(i)) && varmax.vn(i) <= 0
                 varmax.vn(i)=var.v(grd.nx/2,i);
            end
        end
        plot(grd.yv, var.v(grd.nx/2,:), grd.yv, varmax.vn(:), grd.yv, varmax.vp(:))
        title('v')
        xlabel('y (x=0.5)')
        
        subplot(2,2,4)
        var.gy=abs(calc_grady(var.v, grd));
        error=var.gy(grd.nx/2,:); %-var0.ke(grd.nx/2,:);
        for i=1:grd.ny
            if (error(i)>varmax.qp(i)) && varmax.qp(i) >= 0
                varmax.qp(i)=error(i);
            elseif (error(i)<varmax.qn(i)) && varmax.qn(i) <= 0
                 varmax.qn(i)=error(i);
            end
        end
        %plot(grd.yh, error);
        plot(grd.yz, error, grd.yz, ...
            varmax.qn, grd.yz, varmax.qp)
        title('ke')
        xlabel('y (x=0.5)')

%         subplot(2,2,4)
%         error=var.q(grd.nx/2,:); %-var0.q(grd.nx/2,:);
%         for i=1:grd.ny
%             if (error(i)>varmax.qp(i)) && varmax.qp(i) >= 0
%                 varmax.qp(i)=error(i);
%             elseif (error(i)<varmax.qn(i)) && varmax.qn(i) <= 0
%                  varmax.qn(i)=error(i);
%             end
%         end
%         plot(grd.yz, error);
%         %plot(grd.yz, error, grd.yz, ...
%         %    varmax.qn, grd.yz, varmax.qp)
%         title('q')
%         xlabel('y (x=0.5)')
        %title(['Energy var t = ', num2str(t,'%6.2f')])
        %xlabel('time step')
        
%         
%         subplot(3,2,5)
%         for i=1:grd.ny
%             if (var.tmpx(grd.nx/2,i)>varmax.ntxp(i)) && varmax.ntxp(i) >= 0
%                 varmax.ntxp(i)=var.tmpx(grd.nx/2,i);
%             elseif (var.tmpx(grd.nx/2,i)<varmax.ntxn(i)) && varmax.ntxn(i) <= 0
%                  varmax.ntxn(i)=var.tmpx(grd.nx/2,i);
%             end
%         end
%         plot(grd.yu, var.tmpx(grd.nx/2,:), grd.yu, varmax.ntxn(:), grd.yu, varmax.ntxp(:))
%         title('-zv+grad_x(K)')
%         xlabel('y (x=0.5)')
%         
%         subplot(3,2,6)
%         for i=1:grd.ny
%             if (var.tmpy(grd.nx/2,i)>varmax.ntyp(i)) && varmax.ntyp(i) >= 0
%                 varmax.ntyp(i)=var.tmpy(grd.nx/2,i);
%             elseif (var.tmpy(grd.nx/2,i)<varmax.ntyn(i)) && varmax.ntyn(i) <= 0
%                  varmax.ntyn(i)=var.tmpy(grd.nx/2,i);
%             end
%         end
%         plot(grd.yv, var.tmpy(grd.nx/2,:), grd.yv, varmax.ntyn(:), grd.yv, varmax.ntyp(:))
%         title('zu+grad_y(K)')
%         xlabel('y (x=0.5)')
%         
        drawnow
        
end