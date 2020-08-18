function yout = DS_seperate(v,xdata)
x =v;
global meff
%xdata = frequenzrange;
global cocker_DS
if cocker_DS
    N = x(1)*1e23; %in 10^17/(cm*cm*cm)
    tau = x(2)*1e-15; %in fs
    c = x(3);
    tau2 =x(4)*1e-15; %in fs
else
    N = x(1)*1e23; %in 10^17/(cm*cm*cm)
    tau = x(2)*1e-15; %in fs
    c = x(3);
end
omega = xdata*1e12*2*pi;
echarge = 1.6020e-19;
nenner = 1-(1i.*omega.*tau);
zaehler = N.*echarge.*echarge.*tau./meff;
if cocker_DS
    nenner2 = 1-(1i.*omega.*tau2);   
    DS1 = (zaehler./nenner) .*(1 + (c./nenner2));
else
    DS1 = (zaehler./nenner) .*(1 + (c./nenner));
end

yout(:,1) = real(DS1);
yout(:,2) = imag(DS1);