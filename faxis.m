% [f df N dt] = faxis (t, N)
% Berechnet Frequenzachse mit N Punkten aus gegebener Zeitachse.
function [f df Nt dt] = faxis (t, N)
Nt = length (t);
if ~exist ('N', 'var'), N = Nt; end
dt = (t(Nt)-t(1))/(Nt-1);
df = 1/N/dt;

% Neu 28.4.10: Bisher anfällig für Rundungsfehler. Daher:
%f = 0:df:1/dt-df;
f = (0:N-1)*df;

% Neu 30.3.10: Sonderbehandlung für 1-Punkt-Achsen
if Nt == 1
    f = zeros (1, N);
end

s = size (t);
if s(1) > s(2)
    f = f';
end
