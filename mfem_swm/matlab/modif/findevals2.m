% To read in the system matrix in file aaa.m
% and find and plot the eigenvalues
%
% This version assumes the system matrix is for the result
% of one time step rather than the time derivatives

clear

% Read in data
%aaa      % u00 = 200
%aaa0    % standard parameters u00 = 20
aaa00   % u00 = 0
a0 = reshape(a,nmat,nmat);

%aaa1 % rotatn inc by 1%
%a1 = reshape(a,nmat,nmat);

%aaa2 % u00 inc by 1%
%a2 = reshape(a,nmat,nmat);

%da1 = a1 - a0;
%da2 = a2 - a0;

% Arrays of true eigenvalues:
% gravity modes and Rossby modes
dt = 60.0;
rearth = 6371220;
rotatn = 7.29212e-5;
phi00 = 1.0e8;
u00 = 0.0;
rot2 = 2.0*(rotatn + u00/rearth);
om0 = sqrt(phi00)/rearth;
% q = ceil(sqrt(nmat))
q = 17;
ixg = 0;
ixr = 1;
ertru(1) = 0.0;
for n = 1:q
    om = sqrt(n*(n+1)*om0*om0);
    for m = -n:1:n
        ixg = ixg+1;
        egtru(ixg) = om;
        %ixg = ixg+1;
        %egtru(ixg) = -om;
        if (n <q | abs(m) < q-1) % Trick to get the right number of modes
            ixr = ixr + 1;
            ertru(ixr) = m*((u00/rearth) - rot2/(n*(n+1)));
        end
    end
end


% Find evals and right evecs
[v,e] = eig(a0);
% and left evecs
er = diag(e);
[w,e] = eig(a0.');
el = diag(e);

% Convert amplification factors to frequencies
% and damping rates
omegar = - imag(log(er))/dt;
rater =  -real(log(er))/dt;

omegal = - imag(log(el))/dt;
ratel =  -real(log(el))/dt;

% Sort according to frequency
[junk,idx] = sort(omegar);
omegar = omegar(idx);
rater = rater(idx);
er = er(idx);
v = v(:,idx);
[junk,idx] = sort(omegal);
omegal = omegal(idx);
ratel = ratel(idx);
el = el(idx);
w = w(idx,:);

% Set index ranges by hand
g1 = [1:161];
g2 = [482:642];
r1 = [162:481];

% Check L and R evals are thee same
% subplot(2,1,1)
% plot(omegar(r1),'+k')
% hold on
% plot(omegal(r1),'ok')
% hold off
% pause


subplot(2,1,1)

plot(omegal(g2),'+k')
hold on
plot(sort(egtru),'ok')
hold off
title('Hexagonal 642 dof      Gravity modes')
ylabel('Frequency')
axis([0 180 0 2e-2])

subplot(2,1,2)
plot(ratel(g2),'+k')
ylabel('Damping rate')

pause


subplot(3,1,1)

plot(omegal(r1),'+k')
hold on
plot(sort(ertru),'ok')
hold off
title('Hexagonal 642 dof      Rossby modes')
ylabel('Frequency')
axis([0 320 -1e-4 1e-4])

subplot(3,1,2)
plot(ratel(r1),'+k')
ylabel('Damping rate')
rmin = min(ratel(r1));
title(['Minimum = ' num2str(rmin) ])
axis([0 320 -2e-5 4e-5])

subplot(3,1,3)
plot(omegal(r1),'+k')
hold on
plot(sort(ertru),'ok')
hold off
title('Lowest frequency Rossby modes')
ylabel('Frequency')
axis([121 200 -5e-6 5e-6])

pause

% Estimate zonal wavenumber from sensitivity
% of eigenvalue to rotation rate and background wind.
%drot = 0.01*rotatn;
%du00 = 0.01*u00;
%for ix = g2
%    e0 = er(ix);
%    vr = v(:,ix);
%    vl = w(ix,:);
%    de1 = (vl*da1*vr) / (vl*vr);
%    de2 = (vl*da2*vr) / (vl*vr);
%    dom1 = -imag(de1/e0)/dt;
%    dom2 = -imag(de2/e0)/dt;
%    domdrot = dom1/drot;
%    domdu = dom2/du00;
%    m = rearth*domdu   % GW
%    m = rearth*domdu - domdrot; % RW
%    pause
%end




