function [chi] = new__Eon_durch_Eoff_no_model_v0(x); % 
global n_off
global d_on
global d_gesamt
global d_off
global frequenzpunkt
global fschritt
global fmitteln
global n_substrat
global Etarget
global FractionE
global drudesmith
global freq_sno
global f_sio2
global freq_long
global E_on_off_long
global internal_FP
global tapecell
global dteflon

Eon_durch_Eoff = E_on_off_long;
freq = freq_long;
n_air = 1+ 0*Eon_durch_Eoff;


if fmitteln %Averaging over the step size, helps with the noise 
   index_lb = find(frequenzpunkt>(freq+fschritt/2),1,'last');
   index_ub = find(frequenzpunkt>(freq-fschritt/2),1,'last');
   freq_loc = mean(freq(index_lb:index_ub));
   Etarget = mean(Eon_durch_Eoff(index_lb:index_ub));
   n_air_loc =  mean(n_air(index_lb:index_ub));    
   index_lb = find(frequenzpunkt>(freq_sno+fschritt/2),1,'last');
   index_ub = find(frequenzpunkt>(freq_sno-fschritt/2),1,'last');
   n_off_loc = mean(n_off(index_lb:index_ub));
   index_lb = find(frequenzpunkt>(f_sio2+fschritt/2),1,'last');
   index_ub = find(frequenzpunkt>(f_sio2-fschritt/2),1,'last');    
   n_substrat_loc =  mean(n_substrat(index_lb:index_ub));
else
    index = find(frequenzpunkt>(freq_sno),1,'last');    
    freq_loc = freq_sno(index);
    Etarget = Eon_durch_Eoff(index);
    n_off_loc = n_off(index);
    n_air_loc =  n_air(index);
    n_substrat_loc =  n_substrat(index);
end

%% Labbook 112 JN (the formular at p112 assumes Substrate,layer,substrate! p123 is for air!
%% Implemented internal etalons, JN_2_23-4

n_on_loc = x(1) - 1i*x(2);
c0 =physconst('LightSpeed');
freq_loc_SI = freq_loc*1e12;

k0 = 2*pi*freq_loc_SI/c0;


if tapecell
     % index explanation: a = air; t = tape; s = sample (n_off); p = sample(n_on) (photoexcited)
    n_tape =n_substrat_loc; n_on =n_on_loc;
    rho_at = (n_air_loc - n_tape)./(n_air_loc + n_tape);  tau_at = 1+rho_at; rho_ta = -rho_at; tau_ta = 1+rho_ta;
    rho_tp = (n_tape-n_on)./(n_tape+n_on);                tau_tp = 1+rho_tp; rho_pt = -rho_tp; tau_pt = 1+rho_pt;
    rho_ts = (n_tape-n_off_loc)./(n_tape+n_off_loc);      tau_ts = 1+rho_ts; rho_st = -rho_ts; tau_st = 1+rho_st;
    rho_sp = (n_off_loc-n_on)./(n_off_loc+n_on);          tau_sp = 1+rho_sp; rho_ps = -rho_sp; tau_ps = 1+rho_ps;
    if internal_FP
        %case 2, labbook JN2 p54
        I = tau_tp/tau_ts * tau_pt/tau_st;
        II = exp(-1i*k0*n_off_loc*(d_off-d_gesamt));
        III = tau_ps * exp(-1i*k0*n_on*d_on);
        IV = (1+rho_at*rho_ts*exp(-2*1i*k0*dteflon*n_tape)) ./(1+rho_at*rho_tp*exp(-2*1i*k0*dteflon*n_tape));
        V = 1./(1+rho_tp*rho_ps*exp(-2*1i*k0*d_on*n_on));
        VI= 1./(1+rho_ps*rho_st*exp(-2*1i*k0*d_off*n_off_loc));
        VII=   (1+rho_ts*rho_st*exp(-2*1i*k0*d_gesamt*n_off_loc)); %from the reference measurement, that why not 1./...
        FractionE = I*II*III*IV*V*VI*VII;
    else % Case1 labbook JN2 p. 53
        I = tau_tp./tau_ts * tau_pt./tau_st;
        II = exp(-1i*k0*d_on*(n_on-n_off_loc));
        III = (1+rho_at*rho_ts*exp(-2*1i*k0*dteflon*n_tape)) ./(1+rho_at*rho_tp*exp(-2*1i*k0*dteflon*n_tape));
        IV =  (1+rho_ts*rho_st*exp(-2*1i*k0*d_gesamt*n_off_loc)) ./(1+rho_tp*rho_pt*exp(-2*1i*k0*d_on*n_on));
        V =   (1+rho_st*rho_ta*exp(-2*1i*k0*dteflon*n_tape)) ./(1+rho_pt*rho_ta*exp(-2*1i*k0*dteflon*n_tape));
        FractionE = I*II*III*IV*V;
    end
else
    
    rho1_off = (n_substrat_loc - n_off_loc)/(n_substrat_loc + n_off_loc);
    rho2_off = ( n_off_loc - n_air_loc)/( n_off_loc + n_air_loc);
    rho1_on = (n_substrat_loc - n_on_loc)/(n_substrat_loc + n_on_loc);
    rho2_on = ( n_on_loc - n_air_loc)/( n_on_loc + n_air_loc);
    tau1_off = 1+rho1_off;
    tau2_off = 1+rho2_off;
    tau1_on = 1+rho1_on;
    tau2_on = 1+rho2_on;
    
    if internal_FP
        %% see labbook JN2_22-24
        I = n_on_loc/n_off_loc;
        II = (n_off_loc+n_air_loc)./(n_on_loc+n_air_loc);
        III = (2*n_off_loc)./(n_off_loc+n_on_loc);
        IV = exp(1i*k0*d_on*(n_off_loc - n_on_loc));
        rho_n_nstern = (n_off_loc - n_on_loc)/(n_off_loc + n_on_loc);
        rho_nster_0 = (n_on_loc - n_air_loc)/(n_on_loc + n_air_loc);
        Va = 1/(1+ (rho_n_nstern * rho_nster_0 * exp(-1i*2*k0*n_on_loc*d_on))); % Etalon in the photoexcited
        rho_s_n = (n_substrat_loc - n_off_loc)/(n_substrat_loc + n_off_loc);
        Vb = 1/(1+ (rho_s_n * rho_n_nstern * exp(-1i*2*k0*n_off_loc*d_off))); % Etalon in the non-excited between substrat and excited
        rho_n_0 = (n_off_loc - n_air_loc)./(n_off_loc + n_air_loc);
        Vc = 1/(1+ (rho_s_n * rho_n_0 * exp(-1i*2*k0*n_off_loc*d_gesamt)));
        FractionE = I*II*III*IV*Va*Vb/Vc;
    else
        delta_exp = exp(-1i*k0*d_on*(n_on_loc-n_off_loc));
        zaehler = 1+(rho1_off*rho2_off *(exp(-1i*2*k0*n_off_loc*d_on)));
        nenner = 1+(rho1_on*rho2_on *(exp(-1i*2*k0*n_on_loc*d_on)));
        %% See Labbook 123 JN
        FractionE = (n_on_loc/n_off_loc) * ((n_substrat_loc + n_off_loc)/(n_substrat_loc + n_on_loc)) * delta_exp * zaehler/nenner *((n_off_loc+1)/(n_on_loc+1));
    end
end
%% Labbook 107; [1]
chi1 = (log(abs(FractionE)) - log(abs(Etarget)))^2;
chi2 = (angle(FractionE) - angle(Etarget))^2;
chi = chi1+chi2;




%[1]Duvillaret, L. and Garet, F. and Coutaz, J. L.
% A reliable method for extraction of material parameters in terahertz time-domain spectroscopy
%IEEE Journal of Selected Topics in Quantum Electronics 1996
%Vol 2, I3
end
