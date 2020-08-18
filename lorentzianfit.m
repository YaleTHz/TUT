function [chi] = lorentzian(x); % 
global omegarange
global epsilon_i
global epsilon_r

offset = x(1);
omegap = x(2);
omega0 = x(3);
gammar = x(4);
oo =omegarange.^2;

nenner = ((omega0^2-oo)).^2 + (oo).*(gammar^2);
zaehlerreal = (omega0^2-oo)*(omegap^2);
zaehlerimag = omegarange.*omegap.*omegap.*gammar;

epsilon_r_l = offset+ (zaehlerreal./nenner);
epsilon_i_l = zaehlerimag./nenner;
realediff = abs(epsilon_r_l-epsilon_r);
imagdiff = abs(epsilon_i_l-epsilon_i);
chichi = realediff+10*imagdiff;
chi = sum(chichi);
end
