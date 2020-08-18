function [n_out,k_out] = retrieval_of_n_airref_nooptimization(n_in,k_in); % 
%% DO NOT EDIT! 
% this function is used by multiple programs. Copy and rename the function
% prior  to editing! 
%Labbook JN 59-61
% adaption to thin layer, page 90 ff
% Teflon stands for fused quartz (SiO2)
% DLN stands for pilA
%x(1)=2; x(2) =1;
global frequenzpunkt
global freq_air  
global freq_sample
global freq_reference
global phasenorm
global specnorm
global mbranch
global Eair
global Ereference
global Esample
global kteflon 
global nteflon
global dteflon
global dpila
global teflon
global S21teflonphase
global S21DLNphase
global frequenzschritt
global frequenzmitteln
global frequenzrange
global datenfensternfuerfit 
datenfensternfuerfit =1;
global S21target
global S21teflon
%%varliste
% global freqstart freqende fft_air_windowed fft_reference_windowed
% fit B1-C3, d
% implement in optimization toolbox! 
f = frequenzpunkt;
if frequenzmitteln 
    range = frequenzschritt/2000;% /1000 for GHz->THz, /2 for half +/-
    argh = f-range;
    [indexair_lb mini] = find(freq_air>argh,1,'first');
    [indexsample_lb mini] = find(freq_sample>argh,1,'first');
    [indexreference_lb mini] = find(freq_reference>argh,1,'first');
    argh = f+range;
    [indexair_ub mini] = find(freq_air>argh,1,'first');
    [indexsample_ub mini] = find(freq_sample>argh,1,'first');
    [indexreference_ub mini] = find(freq_reference>argh,1,'first'); 
end
[indexair mini] = find(freq_air>f,1,'first');
[indexsample mini] = find(freq_sample>f,1,'first');
[indexreference mini] = find(freq_reference>frequenzpunkt,1,'first');
[indexforn minini]= find(frequenzrange<frequenzpunkt,1,'last');
%%constant
c0 = 2.99792458e8;
wellenl = (c0./f)*1e-6;
k0 = (2*pi)./wellenl;

ff = f*1e12;
wellenlaengeSI = c0/ff;

K = (2*pi)/wellenlaengeSI;
k0 =K;
h = dpila*1e-6;
d = dteflon *1e-6;
teflon=0;
%% refractive index
nSiO = nteflon(indexforn)+1i*kteflon(indexforn);

npilA =n_in +1i*k_in;
global S21DLN

if frequenzmitteln 
        S21target  = mean(S21DLN(indexreference_lb:indexreference_ub));
        S21targetphase  = mean(S21DLNphase(indexreference_lb:indexreference_ub));
else
        S21target = S21DLN(indexsample);
        S21targetphase = S21DLNphase(indexsample);
end

if frequenzmitteln 
        S21target_off  = mean(S21teflon(indexreference_lb:indexreference_ub));
        S21targetphase_off  = mean(S21teflonphase(indexreference_lb:indexreference_ub));
else
        S21target_off = S21teflon(indexsample);
        S21targetphase_off = S21teflonphase(indexsample);
end

global testmich
testmich =indexforn;
n0 =1+1i*0;


%% Orfanidis 158, 165, and some own shit from me.... 



m = -25:1:25;


global testmich
testmich =indexreference;
nn = (-((S21targetphase-S21targetphase_off)+(m*2*pi))/(k0*h))+1; %aphase delay for infty reflections
minimirmich = abs((nn-n_in));
[nutzlos indexen] = min(minimirmich);
n_out = nn(indexen);
mbranch(indexreference) = m(indexen);

kmax = (k_in+2);
kmin = (k_in-2);
if kmin < 0
    kmin =0;
end
kstep = (kmax-kmin)/50000;

krange = kmin:kstep:kmax;

n2 = n_out+1i*krange;
npilA = n2;
tau12 = abs((2*n0)/(n0+nSiO)); %from air into SiO2
P = exp(-1i*k0*(imag(nSiO))*d);

rho1 = abs(((nSiO-npilA)/(nSiO+npilA))); %from sio2 to pila
rho2 = abs(((npilA-n0)/(npilA+n0))); %from pilA to air
tau1 = 1+rho1; tau2= 1+rho2;

Zaehler = tau1*tau2*exp(-1i*k0*npilA*h);
nenner =1+(rho1*rho2*exp(-1i*2*k0*npilA*h));

S21range = (Zaehler./nenner) .*tau12.*P;

[nutzloser indexer] = min(abs(abs(S21target)-S21range));
k_out = krange(indexer);






%T = (Zaehler/nenner)*tau12*P;
%chi1 = (abs(T)-abs(S21target))^2;
%chi2 = ((unwrap(angle((T)))-S21targetphase)^2);
%chi =chi1+chi2;

%chi = abs(T-S21target);
%chi = sqrt((real(T)-real(S21target))^2 + (imag(T)-imag(S21target))^2);



end

