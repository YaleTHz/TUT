function [chi] = retrieval_of_n_airref_thinlayer(x); % 
%% DO NOT EDIT! 
% this function is used by multiple programs. Copy and rename the function
% prior  to editing! 
% Labbook JN 59-61
% adaption to thin layer, page 90 ff
% Teflon stands for fused quartz (SiO2)
% DLN stands for pilA
% DLN and/or pila stands for SnO2 
%x(1)=2; x(2) =1;

%% call global
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
global Esample_on
global kteflon 
global nteflon
global dteflon
global dpila
global teflon
global S21teflonphase
global phasecalculated
global reference_on_air
global lit_val_substrate
global lit_val_substrate_value
global frequenzschritt
global frequenzmitteln
global frequenzrange
global datenfensternfuerfit 
datenfensternfuerfit =1;
global S21target
global mmmmm
global S21_on
global S21_off
global S21_onphase
global S21_offphase
global on_or_off
global S21DLN
global S21prime
global noswitch
global output_fit_opti
%% frequency averaging
f = frequenzpunkt; %single point calculation
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
if reference_on_air
    S21 =  S21DLN;
    S21phase = unwrap(angle(S21DLN));
else
    S21 =  S21prime;
    S21phase = unwrap(angle(S21prime));
end

if frequenzmitteln 
        S21target  = mean(S21(indexreference_lb:indexreference_ub));
        S21targetphase  = mean(S21phase(indexreference_lb:indexreference_ub));
else
        S21target = S21(indexsample);
        S21targetphase = S21phase(indexsample);
end

%% constants and wave vectors
c0 = physconst('LightSpeed'); 
wellenl = (c0./f)*1e-6;
k0 = (2*pi)./wellenl;
ff = f*1e12;
wellenlaengeSI = c0/ff;
K = (2*pi)/wellenlaengeSI;
k0 =K; % there was a time, when I had a strange cgs system is this function ... however by now everything should be SI
h = dpila*1e-6;
d = dteflon *1e-6;
teflon=0;
%% refractive index
nSiO = nteflon(indexforn)-1j*kteflon(indexforn); %loads the previous substrate value from the global variables
if lit_val_substrate
    nSiO = nSiO*0+lit_val_substrate_value;
end
npilA =x(1) -1j*x(2); %builds the complex refractive index from n and k
m =0; %branch variable of retriveal 
global testmich
testmich =indexforn; % This will be the output variable for the 2pi branch m 
n0 =1-1i*0; %Air 
%% Orfanidis 158, 165, and some own shit from me.... 
tau12 = (2*n0)/(n0+nSiO); %from air into SiO2
P = exp(-1i*k0*(nSiO-1)*d).*exp(-1i*k0*(-1)*h);%the second exp function accounts for the fact that air is replaced with SnO2

rho1 = ((nSiO-npilA)/(nSiO+npilA)); %from sio2 to pila
rho2 = ((npilA-n0)/(npilA+n0)); %from pilA to air
tau1 = 1+rho1; tau2= 1+rho2; %please note, I'm inconsistant with the indices .... tau1 means tau01 while tau2 is tau12 ....

if reference_on_air
    Zaehler = tau1*tau2*exp(-1j*k0*(npilA)*h); %nominator
    nenner =1+(rho1*rho2*exp(-1j*2*k0*npilA*h)); %denominator
    T = (Zaehler/nenner)*tau12*P;
else
    tau_10 = (2*nSiO)./(nSiO+n0);
    tau_1n = (2*nSiO)./(nSiO+npilA);
    tau_n0 = (2*npilA)./(npilA+n0);
    taus =  tau_1n*tau_n0/tau_10;
    Pn = exp(-1j*k0*(npilA-n0)*h);
    FP = 1./(1+(rho1*rho2*exp(-1j*2*k0*npilA*h)));
    T = taus*Pn*FP;
end

%% Checking if the right branch is chosen
branchtest = angle(T)-angle(S21target)+m*2*pi; %Ensures that m is the right value
if abs(branchtest)> 1.01*pi
    m =1;
end
branchtest = angle(T)-angle(S21target)+m*2*pi;
if abs(branchtest)> 1.01*pi
    m =-1;
end
branchtest = angle(T)-angle(S21target)+m*2*pi;
if abs(branchtest)> 1.01*pi
    m =-2;
end
branchtest = angle(T)-angle(S21target)+m*2*pi;
if abs(branchtest)> 1.01*pi
    m =2;
end
branchtest = angle(T)-angle(S21target)+m*2*pi;
if abs(branchtest)> 1.01*pi
    m =3;
end
branchtest = angle(T)-angle(S21target)+m*2*pi;
if abs(branchtest)> 1.01*pi
    m =-3;
end
branchtest = angle(T)-angle(S21target)+m*2*pi;
if abs(branchtest)> 1.01*pi
    m =0;
end
%% 
m=0;
chi1 = (angle(T)-(angle(S21target)+m*2*pi))^2; 
chi2 = ((log(abs(T)) - log(abs(S21target)))^2); 
chi = chi1+chi2; % Goal function
%chi=chi1;
%chi = sqrt((real(T)-real(S21target))^2 + (imag(T)-imag(S21target))^2);
mmmmm =m;
phasecalculated = angle(T);
output_fit_opti = T;
end

