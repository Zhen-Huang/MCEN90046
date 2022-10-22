clear
s = 7.5;                    % semi span (m)
c = 2;                      % chord (m)
m = 200;                    % mass per unit area (kg/m^2)

rho = 1.225;                % air density  (kg/m^3)
aw = 2 * pi;                  % 2D lift curve slope
Mtdot = -1.2;               % unsteady torsional damping term
%Mtdot = 0;                 % quasi-steady torsional damping term

EI =  2e+7  ;               % flexural rigidity
GJ =  2e+6  ;               % torsional rigidity
xf = 0.48 * c;              % position of elastic axis
e = (xf - 0.25 * c)/c;      % eccentricity ratio between elastic axis and aerodynamic centre
xcg = 0.5 * c;              % position of the mass axis in a chordwise sense

GJstr = sprintf('%0.5g',GJ); EIstr = sprintf('%0.5g',EI);
astring = ['xf/c = ', num2str(xf/c),' xcg/c = ',num2str(xcg/c),' EI = ',EIstr,...
    ' GJ = ',GJstr,' Mtdot = ', num2str(Mtdot)];
disp(astring)

mt = m * (6 * xcg / c - 2);
ml = 2 * m - mt;
a = ml;
b = (mt - ml) / c;

modes = 2;

% Inertia matrix
A = zeros(modes,modes);
A(1,1) = a * s / 5 * c + b * s * c^2 / 10;
A(1,2) = a * s / 4 * (c^2 / 2 - c * xf) + b * s / 4 * (c^3 / 3 - c^2 / 2 * xf);
A(2,1) = A(1,2);
A(2,2) = a * s / 3 * (c^3 / 3 - c^2 * xf + xf^2 * c) + b * s / 3 * (c^4 / 4 - 2 * c^3 * xf / 3 + xf^2 * c^2 / 2);

% Structural damping matrix
D = zeros(modes, modes);

% Stiffness matrix
E = zeros(modes,modes);
E(1,1) = 4 * EI / s^3;
E(1,2) = 0;
E(2,1) = 0;
E(2,2) = GJ / s;

% Eigenvalue solution, normalisation and sorting for wind off system
[evects, evals1] = eig(A \ E) ;       % eigenvalues / vectors
disp(['wind off nat freqs'])
sfreqs = sqrt(diag(evals1))/(2 * pi)  % natural frequencies
[sfreqs,sort_f] = sort(sfreqs);
sfreqs';
evects = evects(:,sort_f');
max_evect = max(abs(evects));
for jj = 1:modes
    evects(1:modes,jj) = evects(1:modes,jj) / max_evect(jj);
end

% Aerodynamic damping matrix
B = zeros(modes, modes);
% B(1,1) = c * aw * s / 10;
% B(1,2) = 0;
% B(2,1) = - c^2 * e * aw * s / 8;
% B(2,2) = -c^3 * Mtdot * s / 24;

% Aerodynamic stiffness matrix
C = zeros(modes, modes);
C(1,1) = 0;
C(1,2) = c * aw * s / 8;
C(2,1) = 0;
C(2,2) = -c^2 * e * aw * s / 6;

% Loop around velocities for aeroelastic system
f = []; d = []; evalr = []; evali = []; ev = [];
vstart = 1; vinc = 0.1; vend = 180;
icount = 0;

% Determine first order complex eigenvalue solution at each velocity
for v = vstart:vinc:vend
    icount = icount +1;
    vel(icount) = v;
    
    Q = [zeros(modes, modes) eye(modes, modes); -A\(rho*v^2*C + E) -A\(rho*v*B + D)];
    evals = eig(Q);
    er = real(evals); ei = imag(evals);
    
    for ii = 1:2*modes
        wrad(ii,1) = abs(evals(ii));
        zeta(ii,1) = -real(evals(ii)) / wrad(ii);
        if imag(evals(ii)) == 0
            wrad(ii,1) = 0;
        end
        whz(ii,1) = wrad(ii,1) / (2*pi);
    end
    
    %  [ei, eisort] = sort(ei);    % sort evals in order of imag part of eval
    %  [er, ersort] = sort(er);    % real part of evals for plotting
    %  evals = evals(ersort);      % overall evals are still sorted together
    
    [whz, wsort] = sort(whz);      % sort frequencies in order
    zeta = zeta(wsort);
    
    f = [f whz]; d = [d zeta*100];
    
    evalr = [evalr er(wsort)]; evali = [evali ei(wsort)];
    ev = [ev evals(wsort)];
end

plot(vel,d)
