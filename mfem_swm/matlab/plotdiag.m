% Plot time series of basic diagnostics

fname = 'fort.43';

m = load(fname);

% Number of steps per day
perday = 96;

day = m(:,1)/perday;
phierr = m(:,2);
pverr = m(:,3);
mass = m(:,4);
ape = m(:,5);
ke = m(:,6);
pz = m(:,7);


subplot(2,2,1)
dm = (mass - mass(1))/mass(1);
plot(day,dm,'k-')
title('Relative mass change')
xlabel('Day')
ylabel('\Delta M / M_0')

subplot(2,2,2)
semilogy(day,phierr,'k-',day,pverr,'k--')
title('Mass and PV tracer discrepancy')
xlabel('Day')
ylabel('Max (\Delta X) / Max (X)')

subplot(2,2,3)
ae = ape + ke;
plot(day,ape,'k:',day,ke,'k--',day,ae,'k-')
title('APE (dot), KE (dash) and AE')
xlabel('Day')
ylabel('Energy')

subplot(2,2,4)
dz = (pz - pz(1))/pz(1);
de = (ae - ae(1))/ae(1);
plot(day,de,'k-',day,dz,'k--')
title('Relative AE and Z change')
xlabel('Day')
ylabel('\Delta Z / Z_0       \Delta E / E_0')

