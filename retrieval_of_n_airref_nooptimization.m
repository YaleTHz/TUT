function [n_out,k_out] = retrieval_of_n_airref_nooptimization(n_in,k_in); % 
%Labbook JN 59-61
%% some global stuff 
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
global S21target
global dteflon
global dDLN
global teflon
global S21teflonphase
global S21DLNphase
global frequenzschritt
global frequenzmitteln

global datenfensternfuerfit 
%% varliste

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
%%constant
c0 = 2.99792458e8;
wellenl = (c0./f)*1e-6;
k0 = (2*pi)./wellenl;

ff = f*1e12;
wellenlaengeSI = c0/ff;
K = (2*pi)/wellenlaengeSI;


%% refractive index
if teflon
    d2=dteflon;
else
    d2=dDLN;
end
global S21teflon
global S21DLN
if frequenzmitteln 
    if teflon
        S21target = mean(S21teflon(indexreference_lb:indexreference_ub));
        S21targetphase = mean(S21teflonphase(indexreference_lb:indexreference_ub));
    else
        S21target  = mean(S21DLN(indexreference_lb:indexreference_ub));
        S21targetphase  = mean(S21DLNphase(indexreference_lb:indexreference_ub));
    end
else
    if teflon
        S21target = S21teflon(indexreference);
        S21targetphase = S21teflonphase(indexreference);
    else
        S21target = S21DLN(indexsample);
        S21targetphase = S21DLNphase(indexsample);
    end
end
m = -25:1:25;


global testmich
testmich =indexreference;
nn = (-(S21targetphase+(m*2*pi))/(k0*d2))+1; %approximation: no reflections occure in measurement range
minimirmich = abs((nn-n_in));
[nutzlos indexen] = min(minimirmich);
n_out = nn(indexen);
mbranch(indexreference) = m(indexen);
n1=1;

kmax = (k_in+0.5);
kmin = (k_in-0.5);
if kmin < 0
    kmin =0;
end
kstep = (kmax-kmin)/50000;

krange = kmin:kstep:kmax;

n2 = n_out+1i*krange;
T12 = abs((2*n1)./(n1+n2)); % assumes substrate in air 
T21 = abs((2*n2)./(n1+n2));
efungtion = exp(-k0*d2*krange); %% Watch the sign, i =-j ! 

S21range = T12.*T21.*efungtion;
[nutzloser indexer] = min(abs(abs(S21target)-abs(S21range)));
k_out = krange(indexer);




end

