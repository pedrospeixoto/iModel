%% Plot y slice
function herror=plot_evol(t, k, var, var0, herror, grd, fig)

        figure(fig)
        
        subplot(1,2,1)
        wave=1:1:grd.ny/2;
        wave=wave*2*pi/grd.ny;
        s=fft(var.h(grd.nx/2, :));
        plot(wave,abs(s(2:grd.ny/2+1)))
        title(['\eta Spectrum t = ', num2str(t,'%6.2f')])
        xlabel('Wave num')
        ylabel('Coef Amp')
        %[smax, imax]=max(abs(s(2:end)));
        
        subplot(1,2,2)
        time=2:k;
        time=t*time/k;
        herror(k)=max(max(var.h-var0.h));
        semilogy(time, herror(2:k))
        title('\eta Error')
        xlabel('Time')
        ylabel('Max Error')
        %set(gca,'ytick',[0, pi/4, pi/2, 3*pi/4, pi])
        %set(gca,'yticklabel',[' 0';' \pi/4';' \pi/2';' 3 \pi/4';' pi'])
        
        drawnow
        
end