%% Plot y slice
function plot_en_evol(t, k, energy, enstrophy, fig)

        figure(fig)

        subplot(1,2,1)
        time=0:t/k:t-t/k;
        plot(time, energy(1:k))
        title(['Kin Energy (t = ', num2str(t,'%6.2f'), ')'])
        xlabel('Time')
        %ylabel('% Variation')
        ylabel('KE')
        
        subplot(1,2,2)
        semilogy(0:k-1, enstrophy(1:k))
        title(['Max Grad KE var t = ', num2str(t,'%6.2f')])
        xlabel('Time step')
        ylabel('Max |Grad KE| ')
        
        drawnow
        
end