function [chi] = retrieval_of_d_mark3(x); % d, B1 B2 C1 C2
d2=x(1); %fit variabel in micron
%n2constante = x(2);
global sellerie
global debyemodel
global cauchy
global lorentzian
global freq_air
global teflon
global freq_for_d_start
global freq_for_d_end
global datenfensternfuerfit
global S21teflon
global S21DLN
freqstart = freq_for_d_start;
freqende = freq_for_d_end;
%%varliste
% global freqstart freqende fft_air_windowed fft_reference_windowed
% fit B1-C3, d
% implement in optimization toolbox!
freqstartindex = find(freq_air>freqstart,1,'first');
freqendeindex = find(freq_air>freqende,1,'first');
f = freq_air(freqstartindex:freqendeindex);

if debyemodel
    eps_inf =x(2);
    eps_res1 =x(3);
    tau1=x(4)*1E-12;
    eps_res2 =x(5);
    tau2=x(6)*1E-12;
    eps_res3 =x(7);
    tau3=x(8)*1E-12;
    eps_res4 =x(9);
    tau4=x(10)*1E-12;
    complexomega = 2*pi*1i* f*1E12;
    epsilon_deb = eps_inf +(eps_res1./(1+complexomega.*tau1))  +(eps_res2./(1+complexomega*tau2))+(eps_res3./(1+complexomega*tau3))+(eps_res4./(1+complexomega*tau4));
    n2 = sqrt(epsilon_deb);
end
if cauchy
    A =x(2); B =x(3); C =x(4); alpha =x(5); beta= x(6); gamma = x(7);
    c0 = 2.99792458e8;
    wellenl = (c0./f)*1e-6; %wavelength in micrometer
    ll = wellenl.*wellenl;
    llll = ll.*ll;
    n2 = A + B./ll+C./llll + 1i.*(alpha.*exp(12400.*(beta.*(1./wellenl- 1./gamma))));
end
if sellerie
    B1 = x(2);C1 = x(3); B2 = x(4); C2 = x(5); absorber = x(6);
    c0 = 2.99792458e8;
    wellenl = (c0./f)*1e-6; %wavelength in micrometer
    ll = wellenl.*wellenl;
    llll = ll.*ll;
    n2 = sqrt(1+ B1*ll./(ll-C1) + B2*ll./(ll-C1)) + absorber*1i;
end
if lorentzian
    eps_inf = x(2);
    omega_1 = x(3)*2*pi; %all in THz
    Omega_1 = x(4) *2*pi;
    gamma_1 = x(5);
    omega = f*2*pi;
    zaehler = (omega_1.^2) - (omega.^2) - 1i*gamma_1*omega;
    epsilon = eps_inf + (Omega_1.^2)./(zaehler);
    n2 = sqrt(epsilon);
end



if teflon
    S21target = S21teflon(freqstartindex:freqendeindex);
else
    S21target = S21DLN(freqstartindex:freqendeindex);
end



%%constant
c0 = 2.99792458e8;
wellenl = (c0./f)*1e-6;
k0 = (2*pi)./wellenl;

%% refractive index
n1 =1* ones (size (f));
%n2 = (n2constante-0.015*1i)* ones (size (f)); %fit variabel
ll = wellenl.*wellenl;
llll = ll.*ll;
if exist('n2') ==1
else
    n2= x(1) - 1i*x(2);
end
%n2 = sqrt(1+((B1*ll)./(ll - C1))+((B2*ll)./(ll - C2)));
n3 =1* ones (size (f));
rho12 = (n1-n2)./(n1+n2);
rho23 = (n2-n3)./(n2+n3);
S11 = zeros(size(f));
S21 = zeros(size(f));
T_I = zeros(size(f));
R_I = zeros(size(f));



for j = 1:length(f)
    T12(j) = ((2*n1(j))/(n1(j)+n2(j)));
    T21(j) = ((2*n2(j))/(n2(j)+n1(j)));
    P0(j) = exp(1i*k0(j)*d2*n1(j));
    P2(j) = exp(-1i*k0(j)*d2*n2(j));
    R12(j) = rho12(j);
    R21(j) = -rho12(j);
    if datenfensternfuerfit
        S21(j) = P0(j)*T12(j)*P2(j)*T21(j) ;
    else
        S21(j) = P0(j)*T12(j)*P2(j)*T21(j) *...
            (1+ (R21(j)*P2(j)*R12(j)*P2(j)) +(R21(j)*P2(j)*R12(j)*P2(j))^2+(R21(j)*P2(j)*R12(j)*P2(j))^3);
    end
end
% chireal = sum(abs(real(S21target)-real(S21)));
% chiimag = sum(abs(imag(S21target)-imag(S21)));
%
% chiphase = angle(S21target)- angle(S21);
% phasediff = sum(abs(exp(1i*chiphase)));
% absdiff = sum(abs(abs(S21target)-abs(S21target)));
% chi = 4*absdiff + phasediff;
chi =  sum(abs(S21target-S21));
d2
end

