% Set axes for unstructured grid contour plot

if (ptype == 'latlong')
    axis([0,2*pi,-0.5*pi,0.5*pi])
    % axis([3.05,3.2,-0.04,0.04])
elseif (ptype == 'sphere ')
    jtplotgrid
else
    disp('Please set the value of ptype')
end