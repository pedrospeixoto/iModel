% Plot the cubic or hexagonal grid

lwidth = 1.3;

% Primal grid
if (nprmx == 4)
    coords = load('primalgrid_cube.dat');
else
    coords = load('primalgrid_hex.dat');
end
[ne,six] = size(coords);

% initialize diagnostics of min and max gid lengths
dmin = 1.0;
dmax = 0.0;


psymbol = '-';
for e = 1:ne
    x1 = coords(e,1:3)';
    x2 = coords(e,4:6)';
    jtrotplot
    % pause
end


% Dual grid
if (nprmx == 4)
    coords = load('dualgrid_cube.dat');
else
    coords = load('dualgrid_hex.dat');
end
[ne,six] = size(coords);

% initialize diagnostics of min and max grid lengths
dmin = 1.0;
dmax = 0.0;

psymbol = 'k-';
for e = 1:ne
    x1 = coords(e,1:3)';
    x2 = coords(e,4:6)';
    jtrotplot
    % pause
end

axis off
