%% Initialize variables
function [var, par, varmax]=initialize_tc(grd,  tc, mtd, vtc)



% Height
h(1:grd.nx, 1:grd.ny)=0;

% Bottom topography
b=h;

% Velocities
u=h;
v=h;

% Kinetic energy
ke=h;

% Divergence
dv=h;

% Vorticity
z=h;
zu=h;
zv=h;

% PV
q=h;

% Gradient of Bernouli potential
gx=h;
gy=h;


% Test case
switch tc
    
    case 0 % Read data from files
        
        dir='initial_conditions/';
        %filestem='polvani_B_256_';
        filestem='polvani_E_256_';
        ext='.txt';
        
        display('Reading initial data...')
        strcat(dir, filestem, '*', ext)
        h_read=load(strcat(dir, filestem, 'h', ext));
        u_read=load(strcat(dir, filestem, 'u', ext));
        v_read=load(strcat(dir, filestem, 'v', ext));
        data_read=importdata(strcat(dir, filestem, 'data', ext));
        display('Initial data read.')
        
        ndata=data_read.data(4);
        % Check if data is compatible
        if ndata ~= 2*grd.ny
            display
            display('Dimension of data read does not match the grid wanted')
            ndata
            grd.ny
            error('Data read should have twice the amount of elements as the grid')
        end
        
        %Initialize variables
        %loop over x gridpoints
        for ix=1:grd.nx
            %loop over y gridpoints
            for iy=1:grd.ny
                %h points
                h(ix, iy)=h_read(2*ix, 2*iy);
                
                % u points
                u(ix, iy)=u_read(2*ix-1, 2*iy);
                                
                % v points
                v(ix, iy)=u_read(2*ix, 2*iy-1);
                
                %Kinetic energy at h points
                ke(ix,iy)=0.5*(u(ix, iy)*u(ix, iy)+v(ix, iy)*v(ix, iy));
                
                %Bottom topography
                b(ix, iy)=0;
            end
        end
        u0=max(max(abs(u)));
        v0=max(max(abs(v)));
        h0=max(max(abs(h)));
        Ru=data_read.data(1);
        Rv=Ru;
        Fu=data_read.data(2);
        Fv=Fu;
        B=data_read.data(3);
        f0=u0/Ru/grd.lx;
        g=u0^2/Fu^2/h0;
        c=sqrt(g*h0);
        
    case 1 % x direction dominant trigonometric balanced flow
        
        f0=10; %10;   % Coriolis parameter
        g=10;    % Gravity
        u0=1; %1;    % Reference velocity
        v0=0;    % Reference velocity
        h0=10; %10;   % Reference height
        
        w=2*pi; %Oscilation frequency
        
        %Initialize variables
        %loop over x gridpoints
        for ix=1:grd.nx
            %loop over y gridpoints
            for iy=1:grd.ny
                %h points
                h(ix, iy)=h0+(u0/w)*(f0/g)*cos(w*grd.yh(iy));
                
                % u points
                u(ix, iy)=u0*sin(w*grd.yu(iy));
                %exp(-((grd.yu(iy)-ly/2)^2)*100);
                
                % v points
                v(ix, iy)=0;
                
                %Kinetic energy at h points
                ke(ix,iy)=u0^2*sin(w*grd.yh(iy))*sin(w*grd.yh(iy))/2;
                
                %Vorticity terms (at their appropriate positionings)
                z(ix, iy)=-w*u0*cos(w*grd.yz(iy))+f0; % at z points
                zu(ix,iy)=(-w*u0*cos(w*grd.yv(iy))+f0)*u0*sin(w*grd.yv(iy)); % at v points
                zv(ix,iy)=0; % at u points
                q(ix,iy)=z(ix,iy)/(10+(u0/w)*(f0/g)*u0*cos(w*grd.yz(iy)));
                
                % Gradient of bernoulli potential at u, v points
                gx(ix,iy)=0;
                gy(ix,iy)=-f0*u0*sin(w*grd.yv(iy))+u0^2*w*sin(w*grd.yv(iy))*cos(w*grd.yv(iy));
                
                % Divergence at h points
                dv(ix, iy)=0;
                
                %Bottom topography
                b(ix, iy)=0;
            end
        end
        
    case 2 % y direction dominant trigonometric balanced flow
        
        f0=10;   % Coriolis parameter
        g=10;    % Gravity
        u0=0;    % Reference velocity
        v0=1;    % Reference velocity
        h0=10;   % Reference height
        
        %Initialize variables
        %loop over x gridpoints
        for ix=1:grd.nx
            %loop over y gridpoints
            for iy=1:grd.ny
                %h points
                h(ix, iy)=10-(1/(2*pi))*(f0/g)*cos(2*pi*grd.xh(ix));
                
                
                % v points
                v(ix, iy)=sin(2*pi*grd.xv(ix));
                
                % u points
                u(ix, iy)=0;
                
                %Kinetic energy at h points
                ke(ix,iy)=sin(2*pi*grd.xh(ix))*sin(2*pi*grd.xh(ix))/2;
                
                %Vorticity terms (at their appropriate positionings)
                z(ix, iy)=2*pi*cos(2*pi*grd.xz(ix))+f0; % at z points
                zv(ix,iy)=(2*pi*cos(2*pi*grd.xu(ix))+f0)*sin(2*pi*grd.xu(ix)); % at u points
                zu(ix,iy)=0; % at v points
                q(ix,iy)=z(ix,iy)/(10+(1/(2*pi))*(f0/g)*cos(2*pi*grd.xz(ix)));
                
                % Gradient of bernoulli potential at u, v points
                gy(ix,iy)=0;
                gx(ix,iy)=-f0*sin(2*pi*grd.xu(ix))+2*pi*sin(2*pi*grd.xu(ix))*cos(2*pi*grd.xu(ix));
                
                % Divergence at h points
                dv(ix, iy)=0;
                
                %Bottom topography
                b(ix, iy)=0;
            end
        end
        
    case 3 % Damm break - Gaussian
        f0=10; %10;   % Coriolis parameter
        g=10;    % Gravity
        u0=1; %1;    % Reference velocity
        v0=0;    % Reference velocity
        h0=1; %10;   % Reference height
        
        %Initialize variables
        %loop over x gridpoints
        for ix=1:grd.nx
            %loop over y gridpoints
            for iy=1:grd.ny
                %h points
                h(ix, iy)=gaussf(grd.xh(ix), grd.yh(iy), grd, h0);
                %h0+exp(-((grd.xh(ix)-grd.lx/2)^2)*40-((grd.yh(iy)-grd.ly/2)^2)*40);
                
                % v points
                v(ix, iy)=0;
                
                % u points
                u(ix, iy)=0;
                
                %Bottom topography
                b(ix, iy)=0;
            end
            
        end
        
      case 4 % thin layer case with x direction dominant trigonometric balanced flow
        
        f0=0; %10;   % Coriolis parameter
        g=10;    % Gravity
        u0=10; %1;    % Reference velocity
        v0=0;    % Reference velocity
        h0=0.01; %10;   % Reference height
        
        w=2*pi; %Oscilation frequency
        
        %Initialize variables
        %loop over x gridpoints
        for ix=1:grd.nx
            %loop over y gridpoints
            for iy=1:grd.ny
                %h points
                h(ix, iy)=h0; %+(u0/w)*(f0/g)*cos(w*grd.yh(iy));
                
                % u points
                u(ix, iy)= u0*sin(w*grd.yu(iy));
                %exp(-((grd.yu(iy)-ly/2)^2)*100);
                
                % v points
                v(ix, iy)=0;
                
                %Kinetic energy at h points
                ke(ix,iy)=u0^2*sin(w*grd.yh(iy))*sin(w*grd.yh(iy))/2;
                
                %Vorticity terms (at their appropriate positionings)
                z(ix, iy)=-w*u0*cos(w*grd.yz(iy))+f0; % at z points
                zu(ix,iy)=(-w*u0*cos(w*grd.yv(iy))+f0)*u0*sin(w*grd.yv(iy)); % at v points
                zv(ix,iy)=0; % at u points
                q(ix,iy)=z(ix,iy)/(h0);
                
                % Gradient of bernoulli potential at u, v points
                gx(ix,iy)=0;
                gy(ix,iy)=u0^2*w*sin(w*grd.yv(iy))*cos(w*grd.yv(iy));
                
                % Divergence at h points
                dv(ix, iy)=0;
                
                %Bottom topography
                b(ix, iy)=(u0/w)*(f0/g)*cos(w*grd.yh(iy));
            end
        end
        
    case 5 % Holingsworth
        %f0=10; %10;   % Coriolis parameter
        %g=10;    % Gravity
        %u0=50/grd.ny; %1;    % Reference velocity
        %v0=0;    % Reference velocity
        %h0=2.5/grd.ny^2; %10;   % Reference height
        f0=0; %grd.ny; %10;   % Coriolis parameter
        g=10;    % Gravity
        u0=10; %500/128; %1;    % Reference velocity
        v0=0;    % Reference velocity
        h0=0.1; %0/128^2; %10;   % Reference height
        
        w=2*pi; %Oscilation frequency
        
        
        %Initialize variables
        %loop over x gridpoints
        for ix=1:grd.nx
            %loop over y gridpoints
            for iy=1:grd.ny
                %h points
                h(ix, iy)=h0; %+(h0/10000)*(sin(100*w*grd.yh(iy))); %+(u0/w)*(f0/g)*cos(w*grd.yh(iy));
                
                % u points
                u(ix, iy)= u0*sin(w*grd.yu(iy));
                %exp(-((grd.yu(iy)-ly/2)^2)*100);
                
                % v points
                v(ix, iy)=0;
                
                %Kinetic energy at h points
                ke(ix,iy)=u0^2*sin(w*grd.yh(iy))*sin(w*grd.yh(iy))/2;
                
                %Vorticity terms (at their appropriate positionings)
                z(ix, iy)=-w*u0*cos(w*grd.yz(iy))+f0; % at z points
                zu(ix,iy)=(-w*u0*cos(w*grd.yv(iy))+f0)*u0*sin(w*grd.yv(iy)); % at v points
                zv(ix,iy)=0; % at u points
                q(ix,iy)=z(ix,iy)/(h0);
                
                % Gradient of bernoulli potential at u, v points
                gx(ix,iy)=0;
                gy(ix,iy)=u0^2*w*sin(w*grd.yv(iy))*cos(w*grd.yv(iy));
                
                % Divergence at h points
                dv(ix, iy)=0;
                
                %Bottom topography
                b(ix, iy)=(u0/w)*(f0/g)*cos(w*grd.yh(iy));
                
            end
        end
        %h(grd.nx/2, grd.ny/2)=h(grd.nx/2, grd.ny/2)+(h0/10000);

        
      case 6 % Constant velocity with forcing
        f0=10; %10;   % Coriolis parameter
        g=10;    % Gravity
        u0=50/grd.ny; %1;    % Reference velocity
        v0=0;    % Reference velocity
        h0=2.5/grd.ny^2; %10;   % Reference height
        
        w=2*pi; %Oscilation frequency
        
        %Initialize variables
        %loop over x gridpoints
        for ix=1:grd.nx
            %loop over y gridpoints
            for iy=1:grd.ny
                %h points
                h(ix, iy)=h0; %+(h0/1000)*(cos(2*w*grd.yh(iy))^10); %+(u0/w)*(f0/g)*cos(w*grd.yh(iy));
                
                % u points
                u(ix, iy)= u0; %*sin(w*grd.yu(iy));
                %exp(-((grd.yu(iy)-ly/2)^2)*100);
                
                % v points
                v(ix, iy)=0;
                
                q(ix, iy)=f0/h0;
                %Bottom topography
                %b(ix, iy)=-(u0)*(f0/g)*grd.yh(iy);
                
            end
        end
        h(grd.nx/2, grd.ny/2)=h0+(h0/1000);
        q(grd.nx/2, grd.ny/2)=f0/(h0+(h0/1000));
        
    case 7 % Advection test
        f0=10; %10;   % Coriolis parameter
        g=10;    % Gravity
        u0=10; %1;    % Reference velocity
        v0=0;    % Reference velocity
        h0=1; %10;   % Reference height
        
        %Initialize variables
        %loop over x gridpoints
        for ix=1:grd.nx
            %loop over y gridpoints
            for iy=1:grd.ny
                %h points
                h(ix, iy)=gaussf(grd.xh(ix), grd.yh(iy), grd, h0);
                %h0+exp(-((grd.xh(ix)-grd.lx/2)^2)*40-((grd.yh(iy)-grd.ly/2)^2)*40);
                
                % v points
                v(ix, iy)=0;
                
                % u points
                u(ix, iy)=u0;
                
                %Bottom topography
                b(ix, iy)=0;
            end
            
        end
        
    case 8 % x direction dominant trigonometric balanced flow
                %f0=10; %10;   % Coriolis parameter
        %g=10;    % Gravity
        %u0=50/grd.ny; %1;    % Reference velocity
        %v0=0;    % Reference velocity
        %h0=2.5/grd.ny^2; %10;   % Reference height
        f0=0; %10;   % Coriolis parameter
        g=10;    % Gravity
        u0=20; %500/128; %1;    % Reference velocity
        v0=0;    % Reference velocity
        h0=0.01; %0/128^2; %10;   % Reference height
        
        w=2*pi; %Oscilation frequency
        
        
        %Initialize variables
        %loop over x gridpoints
        for ix=1:grd.nx
            %loop over y gridpoints
            for iy=1:grd.ny
                %h points
                h(ix, iy)=h0; %+(h0/10000)*(sin(100*w*grd.yh(iy))); %+(u0/w)*(f0/g)*cos(w*grd.yh(iy));
                
                % u points
                u(ix, iy)= u0*sin(w*grd.yu(iy));
                %exp(-((grd.yu(iy)-ly/2)^2)*100);
                
                % v points
                v(ix, iy)=0; %(u0/100000)*cos(w*grd.yu(iy));
                
                %Kinetic energy at h points
                ke(ix,iy)=u0^2*sin(w*grd.yh(iy))*sin(w*grd.yh(iy))/2;
                
                %Vorticity terms (at their appropriate positionings)
                z(ix, iy)=-w*u0*cos(w*grd.yz(iy))+f0; % at z points
                zu(ix,iy)=(-w*u0*cos(w*grd.yv(iy))+f0)*u0*sin(w*grd.yv(iy)); % at v points
                zv(ix,iy)=0; % at u points
                q(ix,iy)=z(ix,iy)/(h0);
                
                % Gradient of bernoulli potential at u, v points
                gx(ix,iy)=0;
                gy(ix,iy)=u0^2*w*sin(w*grd.yv(iy))*cos(w*grd.yv(iy));
                
                % Divergence at h points
                dv(ix, iy)=0;
                
                %Bottom topography
                b(ix, iy)=(u0/w)*(f0/g)*cos(w*grd.yh(iy));
                
            end
        end
        %h(grd.nx/2, grd.ny/2)=h(grd.nx/2, grd.ny/2)+(h0/10000);
        
        
        case 9 % x direction dominant trigonometric balanced flow - matches sweet scenario 3
        
        f0=10; %10;   % Coriolis parameter
        g=10;    % Gravity
        u0=1; %1;    % Reference velocity
        v0=0;    % Reference velocity
        h0=10; %10;   % Reference height
        
        w=2*pi; %Oscilation frequency
        
        %Initialize variables
        %loop over x gridpoints
        for ix=1:grd.nx
            %loop over y gridpoints
            for iy=1:grd.ny
                %h points
                h(ix, iy)=h0+sin(w*grd.yh(iy));
                
                % u points
                u(ix, iy)=-g*w*cos(w*grd.yu(iy))/f0;
                %exp(-((grd.yu(iy)-ly/2)^2)*100);
                
                % v points
                v(ix, iy)=0;
                
                %Kinetic energy at h points
                ke(ix,iy)=0.0; %u0^2*sin(w*grd.yh(iy))*sin(w*grd.yh(iy))/2;
                
                %Vorticity terms (at their appropriate positionings)
                z(ix, iy)=0; %-w*u0*cos(w*grd.yz(iy))+f0; % at z points
                zu(ix,iy)=0; %(-w*u0*cos(w*grd.yv(iy))+f0)*u0*sin(w*grd.yv(iy)); % at v points
                zv(ix,iy)=0; % at u points
                q(ix,iy)=0; %z(ix,iy)/(10+(u0/w)*(f0/g)*u0*cos(w*grd.yz(iy)));
                
                % Gradient of bernoulli potential at u, v points
                gx(ix,iy)=0;
                gy(ix,iy)=0; %-f0*u0*sin(w*grd.yv(iy))+u0^2*w*sin(w*grd.yv(iy))*cos(w*grd.yv(iy));
                
                % Divergence at h points
                dv(ix, iy)=0;
                
                %Bottom topography
                b(ix, iy)=0;
            end
        end

end

if tc ~= 0
    %Non dimansional parameters
    Ru=2*u0/grd.dy/f0;    % Rossby Number in u direction
    Rv=2*v0/grd.dy/f0;    % Rossby Number in v direction
    c=sqrt(g*h0);       % Gravity wave speed
    Fu=u0/c;            % Froude Number in u direction
    Fv=v0/c;            % Froude Number in v direction
    B=c^2/f0^2/grd.lx^2; %Burgers number
end
cfl=grd.dt*(u0/grd.dx+v0/grd.dy);      % CFL number

par=struct(...
    'tc', tc, ...  
    'mtd', mtd, ...
    'vtc', vtc, ...
    'f0', f0, ... 
    'g', g, ...  
    'u0', u0, ...
    'v0', v0, ...
    'h0', h0, ...
    'Ru', Ru, ...
    'Rv', Rv, ...
    'Fu', Fu, ...
    'Fv', Fv, ...
    'B', B, ...
    'c', c, ...
    'cfl', cfl ...
    );   

var=struct(...
    'u', u, ...
    'v', v, ...
    'h', h, ...
    'b', b, ...
    'z', z, ...
    'q', q, ...
    'dv', dv, ...
    'ke', ke, ...
    'zu', zu, ...
    'zv', zv, ...
    'gx', gx, ...
    'gy', gy, ...
    'tmpx', h, ...
    'tmpy', h ...
    );

hmax(1:grd.ny)=0.;

varmax=struct(...
    'up', hmax, ...
    'vp', hmax, ...
    'hp', hmax, ...
    'qp', hmax, ...
    'ntxp', hmax, ...
    'ntyp', hmax, ...
    'un', hmax, ...
    'vn', hmax, ...
    'hn', hmax, ...
    'qn', hmax, ...
    'ntxn', hmax, ...
    'ntyn', hmax ...
    );


end