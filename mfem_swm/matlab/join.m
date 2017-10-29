% join two points on contour plot

longo1 = long(ifo1);
lato1 = lat(ifo1);
longo2 = long(ifo2);
lato2 = lat(ifo2);

long1 = long(if1);
lat1 = lat(if1);
long2 = long(if2);
lat2 = lat(if2);

option = 1;
%if (strcmp(ptype, 'latlong'))
%    option = 1
%elseif (strcmp(ptype, 'sphere'))
%    option = 2
%end

if (option == 1)
    x = [betao*longo1 + alphao*longo2, beta*long1 + alpha*long2];
    y = [betao*lato1  + alphao*lato2 , beta*lat1  + alpha*lat2 ];
    dl1 = abs(longo2 - longo1);
    dl2 = abs(long2 - long1);
    dl3 = abs(x(2) - x(1));
    if (dl1 < pi & dl2 < pi & dl3 < pi)
      plot(x,y,'k','linewidth',lwidth)
      hold on
    end
elseif (option == 2)
    temp = cos(lato1);
    xo1 = temp*cos(longo1);
    yo1 = temp*sin(longo1);
    zo1 = sin(lato1);
    temp = cos(lato2);
    xo2 = temp*cos(longo2);
    yo2 = temp*sin(longo2);
    zo2 = sin(lato2);
    x1(1) = betao*xo1 + alphao*xo2;
    x1(2) = betao*yo1 + alphao*yo2;
    x1(3) = betao*zo1 + alphao*zo2;
    temp = cos(lat1);
    xn1 = temp*cos(long1);
    yn1 = temp*sin(long1);
    zn1 = sin(lat1);
    temp = cos(lat2);
    xn2 = temp*cos(long2);
    yn2 = temp*sin(long2);
    zn2 = sin(lat2);
    x2(1) = beta*xn1 + alpha*xn2;
    x2(2) = beta*yn1 + alpha*yn2;
    x2(3) = beta*zn1 + alpha*zn2;
    jtrotplot
else
    disp('Please set the variable ptype')
end