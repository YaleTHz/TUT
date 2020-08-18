function [sigma_ref_air,sigma_ref_SiO2] = retrieval_of_sigma(n_SiO2); % 
%% DO NOT EDIT! 
% this function is used by multiple programs. Copy and rename the function
% prior  to editing! 
% Labbook JN 59-61
% adaption to thin layer, page 90 ff
% Teflon stands for fused quartz (SiO2)
% DLN stands for pilA
%x(1)=2; x(2) =1;
%%
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

%% refractive index
if frequenzmitteln 
        Epila_loc = mean(Esample(indexreference_lb:indexreference_ub));
        Eref_loc = mean(Ereference(indexreference_lb:indexreference_ub));
        Eair_loc = mean(Eair(indexreference_lb:indexreference_ub));
        T = (Epila_loc./Eair_loc);
        t = (Epila_loc./Eref_loc);
else
        Epila_loc = Esample(indexsample);
        Eref_loc = Ereference(indexsample);
        Eair_loc = Eair(indexsample);
        T = (Epila_loc./Eair_loc);
        t = (Epila_loc./Eref_loc);
end
%% conj() forms the complex conjugated. This is needed because matlab has a negative phaser, therefore if physics formulars are used, the phase needs to be complex conjugated. 

%% Labbook, JN, p95-97 (2017)
%% [1]Liang et al., Terahertz Characterization of Single-Walled Carbon Nanotube and Graphene On-Substrate Thin Films;IEEE Transactions on Microwave Theory and Techniques 59 10 (2011)
%% [2]Jaeckl et al.,Microwave properties and strain?induced lattice defects of c?axis?oriented YBa2Cu3O7?? thin films on silicon J Appl Phys, 80(6) 15 (1996)
%% [3]Nuss et al., Terahertz surface impedance of thin YBa2Cu3O7 superconducting films, APL, 58 (22) 1991
%% [4]Walther et al.,Terahertz conductivity of thin gold films at the metal-insulator percolation transition, PR B 76 125408 (2007)

Z0 = 376.73;
h = dpila *1e-6;
d = dteflon *1e-6;

P = (exp(-1i*(n_SiO2-1)*k0*d ));
%n_SiO2 = real(n_SiO2);
sigma_ref_SiO2 =(1/(Z0*h)) * ((1+n_SiO2/t) -1 - n_SiO2);%[2,4]
sigma_ref_air = (1/(Z0*h)) * (((4*n_SiO2)/(T*(n_SiO2+1)))*P -1 -n_SiO2) ;%[1]



end

