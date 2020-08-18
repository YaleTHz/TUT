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
% adaption to tape cell, JN2 48-51

%% call global
global frequenzpunkt
global freq_air  
global freq_sample
global freq_reference
global dteflon
global dpila
global teflon
global phasecalculated
global reference_on_air
global lit_val_substrate
global lit_val_substrate_value
global frequenzschritt
global frequenzmitteln
global datenfensternfuerfit 
datenfensternfuerfit =1;
global S21target
global mmmmm
%global S21_on
%global S21_off
%global S21_onphase
%global S21_offphase
%global on_or_off
global S21DLN
global S21prime
%global noswitch
global n_substrat_loc
%% frequency averaging
f = frequenzpunkt; %single point calculation
if frequenzmitteln
    range = frequenzschritt/2000;% /1000 for GHz->THz, /2 for half +/-
    argh = f-range;
   % [indexair_lb mini] = find(freq_air>argh,1,'first');
   % [indexsample_lb mini] = find(freq_sample>argh,1,'first');
    [indexreference_lb,~] = find(freq_reference>argh,1,'first');
    argh = f+range;
    %[indexair_ub mini] = find(freq_air>argh,1,'first');
   % [indexsample_ub mini] = find(freq_sample>argh,1,'first');
    [indexreference_ub,~] = find(freq_reference>argh,1,'first');
end

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
d_sample = dpila*1e-6;
d_tape = (dteflon *1e-6)/2;
teflon=0;
%% refractive index
 %loads the previous substrate value from the global variables
if lit_val_substrate
    n_tape = n_tape*0+lit_val_substrate_value;
else
    n_tape = n_substrat_loc;
end
n_sample =x(1) -1j*x(2); %builds the complex refractive index from n and k
m =0; %branch variable of retriveal 
global testmich
testmich =n_substrat_loc; % This will be the output variable for the 2pi branch m 
n0 =1-1i*0; %Air 
%% Orfanidis 158, 165, Labbook2 page 48-50
%% Fresnel coeeficients
rho_ts = (n_tape - n_sample)/(n_tape + n_sample);
rho_st = -rho_ts;
rho_at = (n0-n_tape)/(n0+n_tape);
rho_ta = -rho_at;
tau_at = 1+rho_at; tau_ts = 1+rho_ts; tau_st = 1+rho_st; tau_ta = 1+rho_ta;
%% Propagation operator
P_tape = exp(-1i*k0*d_tape*(n_tape-1));
P_sample =exp(-1i*k0*d_sample*(n_sample-1));
%% Fabry-Perot
FP_ats = 1/(1+rho_at*rho_ts*exp(-2*1i*k0*d_tape*(n_tape)));
FP_tst = 1/(1+rho_ts*rho_st*exp(-2*1i*k0*d_sample*(n_sample)));
FP_sta = 1/(1+rho_st*rho_ta*exp(-2*1i*k0*d_tape*(n_tape)));

if reference_on_air
    T = tau_at*P_tape*tau_ts*P_sample*tau_st*P_tape*tau_ta * FP_ats*FP_tst*FP_sta;
else
%not yet implemented
msgbox('the method reference_on_air must be activated! Using a reference material instead of air meassurement is not supported for tape cells, because it is a pain to do ... feel free to try it :)')
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
chi2 = ((log(abs(T)) - log(abs(S21target)))^2)*2; 
chi = chi1+chi2; % Goal function
%chi=chi1;
%chi = sqrt((real(T)-real(S21target))^2 + (imag(T)-imag(S21target))^2);
mmmmm =m;
phasecalculated = angle(T);
global output_fit_opti
output_fit_opti =T;

end

