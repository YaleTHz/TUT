function [chi] = retrieval_of_n_airref_tapecell(x); % 
%% DO NOT EDIT! 
% this function is used by multiple programs. Copy and rename the function
% prior  to editing! 
% Labbook JN 59-61
% adaption to thin layer, page 90 ff
% Teflon stands for fused quartz (SiO2)
% DLN stands for pilA
% DLN and/or pila stands for SnO2 
%x(1)=2; x(2) =1;
% prev. version retrieval_of_n_airref_thinlayer
%% call global
global frequenzpunkt 
global freq_sample
global freq_reference
global dteflon
global phasecalculated
global frequenzschritt
global frequenzmitteln
global frequenzrange
global datenfensternfuerfit 
datenfensternfuerfit =1;
global S21target_teflon
global mmmmm
global S21teflon

%% frequency averaging
f = frequenzpunkt; %single point calculation
if frequenzmitteln 
    range = frequenzschritt/2000;% /1000 for GHz->THz, /2 for half +/-
    argh = f-range;
    [indexreference_lb,~] = find(freq_reference>argh,1,'first');
    argh = f+range;
    [indexreference_ub,~] = find(freq_reference>argh,1,'first'); 
end
[indexsample,~] = find(freq_sample>f,1,'first');
[indexforn,~]= find(frequenzrange<frequenzpunkt,1,'last');

S21 =  S21teflon;
S21phase = unwrap(angle(S21teflon));

if frequenzmitteln 
        S21target  = mean(S21(indexreference_lb:indexreference_ub));
        S21targetphase  = mean(S21phase(indexreference_lb:indexreference_ub));
else
        S21target = S21(indexsample);
        S21targetphase = S21phase(indexsample);
end
%% constants and wave vectors
c0 = physconst('LightSpeed');
f = frequenzpunkt;
ff = f*1e12;
wellenlaengeSI = c0/ff;
K = (2*pi)/wellenlaengeSI;
k0 =K; % there was a time, when I had a strange cgs system is this function ... however by now everything should be SI
d = dteflon *1e-6;
%% refractive index
n_tape =x(1) -1j*x(2); %builds the complex refractive index from n and k
m =0; %branch variable of retriveal 
global testmich
testmich =indexforn; % This will be the output variable for the 2pi branch m 
n_0 =1-1i*0; %Air 
%% Orfanidis 158, 165, and some own shit from me.... 
rho_ta = (n_tape-n_0)/(n_tape+n_0);
rho_at =-rho_ta;
tau_ta = 1+rho_ta; 
tau_at = 1+rho_at; 
FP_ata = 1/(1+(rho_at*rho_ta*exp(-2*1i*n_tape*d*k0)));
T = (tau_at*exp(-1i*(n_tape-1)*d*k0) *tau_ta * FP_ata);
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
chi1 = (unwrap(angle(T))-unwrap((angle(S21target)+m*2*pi)))^2; 
chi2 = ((log(abs(T)) - log(abs(S21target)))^2)*2; 
chi = chi1+chi2; % Goal function
%chi=chi1;

mmmmm =m;
phasecalculated = angle(T);
S21target_teflon =T;
end

