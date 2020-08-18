function [chi] = new__Eon_durch_Eoff_DS_v0(x); %
global Eon_durch_Eoff
global n_off
global n_air
global freq
global d_on
global d_gesamt
global d_off
global n_substrat
global meff
global fmin
global fmax
global freq_sno
global fminDS
global fmaxDS
global internal_FP
global fit_ds_transfer
global tapecell
global dteflon
global cocker_DS
index_lb = find(freq_sno>fminDS,1,'first');
index_ub =find(freq_sno>fmaxDS,1,'first');




%% Physical constants, I recoomend: do not change
echarge = 1.602*1e-19;
epsilon_0 = 8.854 *1e-12;
c0 =physconst('LightSpeed');
m_e = 9.1093829*1e-31;

%% Drude parameter  
omega = 2*pi*freq_sno*1e12;
k0 = omega/c0;
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
%% Labbook 112 JN (the formular at p112 assumes Substrate,layer,substrate! p123 is for air!
epsilon_l = n_off.^2; %lattice parameter; maybe times 0+1 ...
dn = 1- 1i*omega*tau;
omegapp = (N.*echarge.*echarge.*tau./meff);
if cocker_DS
    dnn = 1- 1i*omega*tau2;
    sigma = omegapp./(dn) .*(1+ c./(dnn));
else
    sigma = omegapp./(dn) .*(1+ c./(dn));
end
epsilon_on = epsilon_l + ((1i*sigma)./(epsilon_0*omega)); %the sign here is linked to the conj in the next line
n_on = sqrt((conj(epsilon_on))); %complex conjugated, inconsitend sign convention.

if tapecell
    % index explanation: a = air; t = tape; s = sample (n_off); p = sample(n_on) (photoexcited)
    n_tape =  n_substrat;
    rho_at = (n_air - n_tape)./(n_air + n_tape);  tau_at = 1+rho_at; rho_ta = -rho_at; tau_ta = 1+rho_ta;
    rho_tp = (n_tape-n_on)./(n_tape+n_on);        tau_tp = 1+rho_tp; rho_pt = -rho_tp; tau_pt = 1+rho_pt;
    rho_ts = (n_tape-n_off)./(n_tape+n_off);      tau_ts = 1+rho_ts; rho_st = -rho_ts; tau_st = 1+rho_st;
    rho_sp = (n_off-n_on)./(n_off+n_on);          tau_sp = 1+rho_sp; rho_ps = -rho_sp; tau_ps = 1+rho_ps;
    if internal_FP
        %case 2, labbook JN2 p54
        I = (tau_tp./tau_ts) .* (tau_pt./tau_st);
        II = exp(-1i.*k0.*n_off.*(d_off-d_gesamt));
        III = tau_ps .* exp(-1i.*k0.*n_on.*d_on);
        IV = (1+rho_at.*rho_ts.*exp(-2.*1i.*k0.*dteflon.*n_tape)) ./(1+rho_at.*rho_tp.*exp(-2.*1i.*k0.*dteflon.*n_tape));
        V = 1./(1+rho_tp.*rho_ps.*exp(-2.*1i.*k0.*d_on.*n_on));
        VI= 1./(1+rho_ps.*rho_st.*exp(-2.*1i.*k0.*d_off.*n_off));
        VII=   (1+rho_ts.*rho_st.*exp(-2.*1i.*k0.*d_gesamt.*n_off)); %from the reference measurement, that why not 1./...
        FractionE = I.*II.*III.*IV.*V.*VI.*VII;
    else % Case1 labbook JN2 p. 53
        I = tau_tp./tau_ts .* tau_pt./tau_st;
        II = exp(-1i.*k0.*d_on.*(n_on-n_off));
        III = (1+rho_at.*rho_ts.*exp(-2.*1i.*k0.*dteflon.*n_tape)) ./(1+rho_at.*rho_tp.*exp(-2.*1i.*k0.*dteflon.*n_tape));
        IV =  (1+rho_ts.*rho_st.*exp(-2.*1i.*k0.*d_gesamt.*n_off)) ./(1+rho_tp.*rho_pt.*exp(-2.*1i.*k0.*d_on.*n_on));
        V =   (1+rho_st.*rho_ta.*exp(-2.*1i.*k0.*dteflon.*n_tape)) ./(1+rho_pt.*rho_ta.*exp(-2.*1i.*k0.*dteflon.*n_tape));
        FractionE = I.*II.*III.*IV.*V;
    end
else
    rho1_off = (n_substrat - n_off)./(n_substrat + n_off);
    rho2_off = ( n_off - n_air)./( n_off + n_air);
    rho1_on = (n_substrat - n_on)./(n_substrat + n_on);
    rho2_on = ( n_on - n_air)./( n_on + n_air);
    tau1_off = 1+rho1_off;
    tau2_off = 1+rho2_off;
    tau1_on = 1+rho1_on;
    tau2_on = 1+rho2_on;
    if internal_FP
        %% see labbook JN2_22-24
        I = n_on./n_off;
        II = (n_off+n_air)./(n_on+n_air);
        III = (2.*n_off)./(n_off+n_on);
        IV = exp(1i.*k0.*d_on.*(n_off - n_on));
        rho_n_nstern = (n_off - n_on)./(n_off + n_on);
        rho_nster_0 = (n_on - n_air)./(n_on + n_air);
        Va = 1./(1+ (rho_n_nstern .* rho_nster_0 .* exp(-1i.*2.*k0.*n_on.*d_on))); % Etalon in the photoexcited
        rho_s_n = (n_substrat - n_off)./(n_substrat + n_off);
        Vb = 1./(1+ (rho_s_n .* rho_n_nstern .* exp(-1i.*2.*k0.*n_off*d_off))); % Etalon in the non-excited between substrat and excited
        rho_n_0 = (n_off - n_air)./(n_off + n_air);
        Vc = 1./(1+ (rho_s_n .* rho_n_0 .* exp(-1i.*2.*k0.*n_off.*d_gesamt)));
        FractionE = I.*II.*III.*IV.*Va.*Vb./Vc;
    else
        delta_exp = exp(-1i.*k0.*d_on.*(n_on-n_off));
        zaehler = 1+(rho1_off.*rho2_off .*(exp(-1i.*2.*k0.*n_off.*d_on)));
        nenner = 1+(rho1_on.*rho2_on .*(exp(-1i.*2.*k0.*n_on.*d_on)));
        %% See Labbook 123 JN
        FractionE = (n_on./n_off) .* ((n_substrat + n_off)./(n_substrat + n_on)) .* delta_exp .* zaehler./nenner .*((n_off+1)./(n_on+1));
    end
end


%% Goal Function, Labbook 107; [1]
% this is the goal function I used for without model
chi1 = sum(sum((log(abs(FractionE(index_lb:index_ub))) - log(abs(Eon_durch_Eoff(index_lb:index_ub)))).^2)); %limiting the fitting range on fmin to max
chi2 = sum(sum((angle(FractionE(index_lb:index_ub)) - angle(Eon_durch_Eoff(index_lb:index_ub))).^2));
chi = chi1+chi2;

% this is the new goal function
%chi1 = sum(sum((real(FractionE(index_lb:index_ub))- real(Eon_durch_Eoff(index_lb:index_ub))).^2));
%chi2 = sum(sum((imag(FractionE(index_lb:index_ub))- imag(Eon_durch_Eoff(index_lb:index_ub))).^2));
%chi = chi1+chi2;

%[1]Duvillaret, L. and Garet, F. and Coutaz, J. L.
% A reliable method for extraction of material parameters in terahertz time-domain spectroscopy
%IEEE Journal of Selected Topics in Quantum Electronics 1996
%Vol 2, I3
fit_ds_transfer = FractionE;
end
