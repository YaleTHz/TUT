function [nDLN_effectivmedium,kDLN_effectivmedium,nDLN_referenceT,kDLN_referenceT,nDLN_effectivemediumT,kDLN_effectivemediumT,...
    nDLN_effectivmediumB,kDLN_effectivmediumB,nDLN_effectivemediumBT,kDLN_effectivemediumBT] =...
    retrieval_of_n_relativtoteflon_nooptimization(nteflon,kteflon,nDLN,kDLN,nDLN_referenceT_old,kDLN_referenceT_old) % 
%Labbook JN 59-61
global frequenzpunkt
global freq_reference
global dteflon
global dDLN
global frequenzschritt
global frequenzmitteln
%%varliste
% global freqstart freqende fft_air_windowed fft_reference_windowed
% fit B1-C3, d
% implement in optimization toolbox! 
f = frequenzpunkt;
if frequenzmitteln 
    range = frequenzschritt/2000;% /1000 for GHz->THz, /2 for half +/-
    argh = f-range;
    [indexreference_lb,~] = find(freq_reference>argh,1,'first');
    argh = f+range;
    [indexreference_ub,~] = find(freq_reference>argh,1,'first'); 
else
[indexreference mini] = find(freq_reference>f,1,'first');
end
%%constant
c0 = 2.99792458e8;
wellenl = (c0./f)*1e-6;
k0 = (2*pi)./wellenl;

%% refractive index
global S21DLN_on_Teflon
global S21DLN_on_Teflon_phase


if frequenzmitteln 
        S21target  = mean(S21DLN_on_Teflon(indexreference_lb:indexreference_ub));
        S21targetphase  = mean(S21DLN_on_Teflon_phase(indexreference_lb:indexreference_ub));
else
        S21target = S21DLN_on_Teflon(indexreference);
        S21targetphase = S21DLN_on_Teflon_phase(indexreference);
end


%% effektives medium, values referenced on air
%nDLN is the refractive index of the teflonsample with DLN normalized on
%air! nteflon if teflon measurement normalized on air


nDLN_compl = nDLN+1i*kDLN;
nteflon_compl = nteflon+1i*kteflon;
epsilonh = nteflon_compl^2; %host materials
epsilonmean = nDLN_compl^2;% host and dopand (amino acid)  ruppin, 2000 optics comunication 182 273-279
global volumenfaktor

zaehler = volumenfaktor*epsilonmean*epsilonh + 2*volumenfaktor*epsilonh*epsilonh -...
            2*epsilonh*epsilonh + 2*epsilonmean*epsilonh;
nenner = volumenfaktor*epsilonmean - epsilonmean + epsilonh + 2*volumenfaktor*epsilonh;

epsilonMG = zaehler/nenner;
nMG = sqrt(epsilonMG);
nDLN_effectivmedium = real(nMG);
kDLN_effectivmedium = imag(nMG);

eta1= volumenfaktor;
eta2 = 1-volumenfaktor;
epsilon2 = epsilonh;
epsilonm = epsilonmean;
GG = (epsilon2-epsilonm)./(epsilon2+2*epsilonm);
zaehler = (eta1*epsilonm)-(2*eta2.*GG.*epsilonm);
nenner = eta1+(eta2*GG);

epsilonBRUGGE = zaehler./nenner;
nBRUGGE = sqrt(epsilonBRUGGE);
nDLN_effectivmediumB = real(nBRUGGE);
kDLN_effectivmediumB = imag(nBRUGGE);



%% seems to work

%now referenced on teflon:



m = -25:1:25;


n1=1;
nn = (dDLN-dteflon+dteflon*nteflon)/dDLN - (S21targetphase-2*pi*m)/(k0*dDLN);
minimirmich = abs((nn-nDLN_referenceT_old));
[~,indexen] = min(minimirmich);
nmixture = nn(indexen);
n1=1;
%%real part of mixture calculated
%%initializing complex part of mixture
k_in = kDLN_referenceT_old;
kmax = (k_in+0.5);
kmin = (k_in-0.5);
if kmin < 0
    kmin =0;
end
kstep = (kmax-kmin)/50000;

krange = kmin:kstep:kmax;

% something is strange! Fix it!

nm = nmixture+1i*krange;
nteflon_compl = nteflon+1i*kteflon;
T12_m = abs((2*n1)./(n1+nm));
T21_m = abs((2*nm)./(n1+nm));
T12_T = abs((2*n1)./(n1+nteflon_compl));
T21_T=  abs((2*nteflon_compl)./(n1+nteflon_compl));
exponent = -k0 *(krange*dDLN - kteflon*dteflon);
efungtion = exp(exponent);
TTTT = T12_m.*T21_m./(T12_T.*T21_T);
S21range = TTTT.*efungtion;
[~,indexer] = min(abs(abs(S21target)-S21range));
kmixture = krange(indexer);

n_mixture_compl = nmixture +1i*kmixture;
nDLN_referenceT = nmixture;
kDLN_referenceT = kmixture;

%% effektives medium, values referenced on teflon
%nDLN is the refractive index of the teflonsample with DLN normalized on
%air! nteflon if teflon measurement normalized on air

clear nMG
clear nBRUGGE
clear indexer
nDLN_compl = n_mixture_compl;
nteflon_compl = nteflon+1i*kteflon;
epsilonh = nteflon_compl^2; %host materials
epsilonmean = nDLN_compl^2;% host and dopand (amino acid)  ruppin, 2000 optics comunication 182 273-279
zaehler = volumenfaktor*epsilonmean*epsilonh + 2*volumenfaktor*epsilonh*epsilonh -...
            2*epsilonh*epsilonh + 2*epsilonmean*epsilonh;
nenner = volumenfaktor*epsilonmean - epsilonmean + epsilonh + 2*volumenfaktor*epsilonh;

epsilonMG = zaehler/nenner;
nMG = sqrt(epsilonMG);
nDLN_effectivemediumT = real(nMG);
kDLN_effectivemediumT = imag(nMG);





eta1= volumenfaktor;
eta2 = 1-volumenfaktor;
epsilon2 = epsilonh;
epsilonm = epsilonmean;
GG = (epsilon2-epsilonm)./(epsilon2+2*epsilonm);
zaehler = (eta1*epsilonm)-(2*eta2.*GG.*epsilonm);
nenner = eta1+(eta2*GG);

epsilonBRUGGE = zaehler./nenner;
nBRUGGE = sqrt(epsilonBRUGGE);
nDLN_effectivemediumBT = real(nBRUGGE);
kDLN_effectivemediumBT = imag(nBRUGGE);








end

