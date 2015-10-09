%% Gráficos
%% Call as >>  PanelFig(n, x - grid levels, matrix-columns are vars -max erros, ymean, 0.1, 0.1, ['a'; 'b'; 'c'])
function PanelFig(n, x, y1, y2, y3, y4, ref1, ref2, names)
marcadores=['-*r';'-og';'-sb';'-^c';'-+m';'-dy';'-vk';'-pk'];
marcvec=cellstr(marcadores);

%Reference lines

ref(1,1)=ref1;
ref(1,2)=ref2;
for j=2:n
    ref(j, 1)=ref(j-1, 1)/2;
    ref(j, 2)=ref(j-1, 2)/4;
end
        fig=figure(1);
        clf;
       title('Interpolation from cell to edges')
        subplot(2,2,1)
        hold off
        
        semilogy(x(1:n),y1(1:n,1),char(marcvec(1)),'LineWidth',1.3)
        hold on
        semilogy(x(1:n),y1(1:n,2),char(marcvec(2)),'LineWidth',1.3)
        semilogy(x(1:n),y1(1:n,3),char(marcvec(3)),'LineWidth',1.3)
        semilogy(x(1:n),y1(1:n,4),char(marcvec(4)),'LineWidth',1.3)
        semilogy(x(1:n),y1(1:n,5),char(marcvec(5)),'LineWidth',1.3)
        semilogy(x(1:n),y1(1:n,6),char(marcvec(6)),'LineWidth',1.3)
        semilogy(x(1:n),ref(1:n,1),'-', 'Color', [192/256 192/256 192/256],'LineWidth',1.3, 'DisplayName', 'O1')
        semilogy(x(1:n),ref(1:n,2),'-', 'Color', [192/256 192/256 192/256],'LineWidth',1.3, 'DisplayName', 'O2')
        ylabel('Max Error','fontsize',12)
        %legend(names(:),8, 'Location', 'Best');
        
        title('Interpolation from cell to edges - Max Error','fontsize',12 )
                
        subplot(2,2,2)
        hold off
        semilogy(x(1:n),y2(1:n,1),char(marcvec(1)),'LineWidth',1.3)
        hold on
        semilogy(x(1:n),y2(1:n,2),char(marcvec(2)),'LineWidth',1.3)
        semilogy(x(1:n),y2(1:n,3),char(marcvec(3)),'LineWidth',1.3)
        semilogy(x(1:n),y2(1:n,4),char(marcvec(4)),'LineWidth',1.3)
        semilogy(x(1:n),y2(1:n,5),char(marcvec(5)),'LineWidth',1.3)
        semilogy(x(1:n),y2(1:n,6),char(marcvec(6)),'LineWidth',1.3)
        semilogy(x(1:n),ref(1:n,1),'-', 'Color', [192/256 192/256 192/256],'LineWidth',1.3, 'DisplayName', 'O1')
        semilogy(x(1:n),ref(1:n,2),'-', 'Color', [192/256 192/256 192/256],'LineWidth',1.3, 'DisplayName', 'O2')
        ylabel('L2 Error','fontsize',12)
        legend(names(:),8, 'Location', 'BestOutside');
        
        title('Intepolation','fontsize',12 )

        hold off
        
                
end