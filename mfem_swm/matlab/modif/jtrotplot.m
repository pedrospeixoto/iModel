% Check to see if this is max or min grid length
% Rotate view point and plot line segment

        x1rot = rot*x1;
        x2rot = rot*x2;
        if (x1rot(1) > 0 | x2rot(1) > 0)
            lx = [x1rot(2),x2rot(2)]';
            ly = [x1rot(3),x2rot(3)]';
            plot(lx,ly,psymbol,'linewidth',lwidth)
            hold on
        end
