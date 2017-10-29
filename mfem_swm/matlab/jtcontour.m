% Draw a contour plot on an unstructured grid

qmin = min(q);
qmax = max(q);

% Set contour values
% v = linspace(45000,55000,11);
% v = linspace(0,1000,11);
% v = linspace(-200.0,200.0,11);
% v = linspace(5000,6000,21);
v = linspace(-3e-9,3e-9,31);
% v = linspace(qmin,qmax,11);
% lw = [2,1,1,1,1,2,1,1,1,1,2];
% lw = [2,2,2,2,2,2,1,1,1,1,1];
lw = [1,1,1,1,1,2,1,1,1,1,1,1,1,1,1,2,...
        1,1,1,1,1,1,1,1,1,2,1,1,1,1,1];

nv = numel(v);


% Find useful range of contours
cstart = 1;
cend = nv;
for iv = 1:nv
    if (v(iv) < qmin)
        cstart = iv;
    end
end
for iv = nv:-1:1
    if (v(iv) > qmax)
        cend = iv;
    end
end


% Set viewing angle for spherical view
alpha = 1.0*pi;     % Longitude of viewpoint
beta = 0.0*pi;         % Latitude of viewpoint
ca = cos(alpha);
sa = sin(alpha);
cb = cos(beta);
sb = sin(beta);
rota = [ca, sa, 0; -sa, ca, 0; 0, 0, 1];
rotb = [cb, 0, sb; 0, 1, 0; -sb, 0, cb];
rot = rota*rotb;

psymbol = 'k-';

% Just to get x1 and x2 the right shape
x1 = [0,0,0]';
x2 = x1;

subplot(1,1,1)

% Loop over polygonal regions
for ir = 1:nr
    
    % Loop over contour values
    for iv = cstart:cend
        
        cv = v(iv);
        lwidth = lw(iv);
        
        % Find edges of region crossed by this contour
        odd = 0;
        hit = 0;
        npr = sum(rlist(ir,:) > 0);
        for ie = 1:npr
            
            iep = ie + 1;
            if (iep > npr)
                iep = 1;
            end
            
            if1 = rlist(ir,ie);
            if2 = rlist(ir,iep);
            q1 = q(if1);
            q2 = q(if2);
                        
            % Does the contour cross this edge?
            d1 = q1 - cv;
            d2 = q2 - cv;
            qq = d1*d2;
            crossed = (qq < 0.0);
            hit = (qq == 0.0 & (d1 + d2 < 0.0));
            if (crossed == 1 | hit == 1)
                                
                odd = 1 - odd;
                if (q1 == q2)
                    alpha = 0.0;
                else
                    alpha = (cv - q1)/(q2 - q1);
                end
                beta = 1.0 - alpha;
                
                % If this is an odd numbered crossing then
                % save it till we find its partner.
                % Otherwise, plot the line
                if (odd == 1)
                    ifo1 = if1;
                    ifo2 = if2;
                    alphao = alpha;
                    betao = beta;
                else
                    join
                end

            end
        end
        if (odd == 1)
            disp(['Region' num2str(ir) 'odd no. of crossings'])
            pause
        end
        
    end
    
end
jtaxes
range = ['  Min ' num2str(min(q)) '  Max ' num2str(max(q))];
title([ytitle  range])
hold off

pause