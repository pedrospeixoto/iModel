function h=gaussf(x, y, grd, h0)
                %h points
                h=h0+exp(-((x-grd.lx/2)^2)*40-((y-grd.ly/2)^2)*40);

end