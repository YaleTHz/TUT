function [chi] = DSfit(x); % 
global n_off
global meff
global fminDS
global fmaxDS
global x_fit
global y_fit
global sigmaa
global cocker_DS

%index_lb = find(x_fit>fminDS,1,'first');
%index_ub =find(x_fit>fmaxDS,1,'first');
index_lb =1;  
%% Physical constants, I recoomend: do not change 
echarge = 1.602*1e-19; 
epsilon_0 = 8.854 *1e-12;
c0 =physconst('LightSpeed');
m_e = 9.1093829*1e-31;

%% Drude parameter
omega = 2*pi*x_fit*1e12;
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

epsilon_l = n_off.^2; %lattice parameter; maybe times 0+1 ...
dn = 1- 1i*omega*tau;
omegapp = (N.*echarge.*echarge.*tau./meff);
if cocker_DS
    dnn = 1- 1i*omega*tau2;
    sigmaa = omegapp./(dn) .*(1+ c./(dnn));
else
    sigmaa = omegapp./(dn) .*(1+ c./(dn));
end
testvalue = sum(isnan(sigmaa));
if testvalue >1
    chi = 1e20;
    msgbox('Error in DSfit.m subroutine')
else
    chi1 = abs((real(y_fit(index_lb:end))-real(sigmaa(index_lb:end))).^2);
    chi2 = abs((imag(y_fit(index_lb:end))-imag(sigmaa(index_lb:end))).^2);
    chi3 = sum(sum(chi1));
    chi4 = sum(sum(chi2));
    chi = abs(chi3+chi4);
end

end