% 3dof_Flutter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plots out the vg and vw plots for a rectangular main wing and full span
% control - 3DOF -  bending / torsion shapes with
% modified strip theory aerodynamics - quasi steady plus unsteady terms
% z = qb(y/s)^2 + qt(y/s)(x-xf) + [x - xh]beta
% Control stiffnesses in book example are kbeta 1e3, 1e4 and 1e5 (Nm/m)/rad
%  IAAL 2nd edition pp193 - 199

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format compact; clear all; close all

% Geometric data
s = 7.5;                    % semi span
c = 2;                      % chord
hline = 80;                 % hinge line as percentage of chord
xh = hline / 100 * c;       % distance of hinge line aft of leading edge
fa = 40;                    % elastic axis as percentage of chord
xf = fa / 100 * c;          % position of elastic axis aft of le
xac = 0.25 * c;             % distance of ac aft of le
xe = xf - xac;              % distance of elastic axis aft of ac
e = xe / c;                 % eccentricity of fa aft of ac per chord

% Mass data
m = 400;                    % mass per unit area
mmain = m;                  % mass per unit area of main surface
mcont = m;                  % mass per unit area of control surface

% Stiffness data
EI =  4e+7;                 % flexural rigidity
GJ =  8e+6;                 % torsional rigidity
kbeta = 1e+4;               % control stiffness per unit span

GJstr = sprintf('%0.5g',GJ); EIstr = sprintf('%0.5g',EI);
astring = ['xf/c = ', num2str(xf/c),' xh/c = ',num2str(xh/c),' EI = ',EIstr,...
    ' GJ = ',GJstr,' kbeta = ', num2str(kbeta)];
disp(astring)

modes = 3;

% Inertia matrix
A = zeros(modes, modes);
% Main surface
Amain = zeros(modes, modes);
Amain(1,1) = xh * s / 5;
Amain(1,2) = (xh^2 / 2 - xh * xf) * s / 4;
Amain(2,1) = Amain(1,2);
Amain(2,2) = (xh^3 / 3 - xh^2 * xf + xf^2 * xh) * s / 3;
% Control surface
Acont = zeros(modes, modes);
Acont(1,1) = (c - xh) * s / 5;
Acont(1,2) = ((c^2 - xh^2) / 2 - xf * (c - xh)) * s / 4;
Acont(1,3) = ((c^2 - xh^2) / 2 - xh * (c - xh)) * s / 3;
Acont(2,2) = ((c^3 - xh^3) / 3 - xf * (c^2 - xh^2) + xf^2 * (c - xh)) * s / 3;
Acont(2,3) = ((c^3 - xh^3) / 3 - xf * (c^2 - xh^2) / 2 - xh * (c^2 - xh^2) / 2 + xf * xh * (c - xh)) * s / 2;
Acont(3,3) = ((c^3 - xh^3) / 3 - xh * (c^2 - xh^2) + xh^2 * (c - xh)) * s;
Acont(2,1) = Acont(1,2);
Acont(3,1) = Acont(1,3);
Acont(3,2) = Acont(2,3);

A = mmain * Amain + mcont * Acont;

% Structural damping matrix is zero
D = zeros(modes, modes);

% Structural stiffness matrix
E = zeros(modes, modes);
E(1,1) = 4 * EI / s^3;
E(2,2) = GJ / s;
E(3,3) = kbeta * s;

% Natural frequencies and mode shapes of the wind off system
[evects, evals1] = eig(A \ E) ;
sfreqs = sqrt(diag(evals1))/(2 * pi);
[sfreqs, sort_f] = sort(sfreqs);
freqs = sfreqs'

% Transformation of normal modes to tip displacements
Ttip = [1^2 -1 * xf 0;
    1^2 1 * (xh - xf) 0;
    1^2 1 * (c - xf) (c - xh)];
Tvects = Ttip * evects;
Tvects = Tvects(:, sort_f');
max_Tvect = max(abs(Tvects));
for jj = 1:modes
    Tvects(1:modes, jj) = Tvects(1:modes, jj) / max_Tvect(jj);
end
Tvects

% Plot mode shapes at wing tip
zline = [0 0 0];
zvals = [0 xh/c 1];
figure(1)
subplot(311)
plot(zvals,Tvects(:,1),'O-')
title([num2str(sfreqs(1),3),'  Hz'])
hold on
plot(zvals,zline,'--b')
subplot(312)
plot(zvals,Tvects(:,2),'O-')
hold on
plot(zvals,zline,'--b')
title([num2str(sfreqs(2),3),'  Hz'])
subplot(313)
plot(zvals,Tvects(:,3),'O-')
hold on
plot(zvals,zline,'--b')
title([num2str(sfreqs(3),3),'  Hz'])

% Bending / torsion / control alone natural frequencies

fbend = 1/2/pi * sqrt(E(1,1)/A(1,1));
ftors = 1/2/pi * sqrt(E(2,2)/A(2,2));
fcont = 1/2/pi * sqrt(E(3,3)/A(3,3));

% Aerodynamic data

% Theodorsen's T functions

chat = 2 * (xh / c) - 1;

T4 = - acos(chat) + chat * sqrt(1 - chat ^ 2);
T10 = sqrt(1 - chat ^ 2) + acos(chat);
T11 = acos(chat) * (1 - 2 * chat) + sqrt(1 - chat ^ 2) * (2 - chat);
T12 = sqrt(1 - chat ^ 2) * (2 + chat) - acos(chat) * (2 * chat + 1);

% Aerodynamic coefficients

aw = 2 * pi;                % lift per incidence
ac = aw * T10 / pi;         % lift per control rotation
bw = e * aw;                % pitching moment per control rotation
bc = e * aw * T10 / pi;     % pitching moment per control rotation
cw = - T12 / 2;             % hinge moment per incidence
cc = - T12 * T10 / 2 / pi;  % hinge moment per control rotation

Mtdot = -1.2;               % unsteady torsional damping term
Mbdot = -0.1;               % unsteady control rotation damping term

rho = 1.225;                % air density
vstart = 1;                 % velocity range
vinc = 1;
vend = 200;

% Aerodynamic matrices

B = zeros(modes, modes);     % aero damping - based on rho*V*B
C = zeros(modes, modes);     % aero stiffness - based on rho*V^2*C

B(1,1) = aw * c * s / 10;
B(1,2) = 0;
B(1,3) = 0;
B(2,1) = - bw * c^2 * s / 8;
B(2,2) = - Mtdot * c^3 * s / 24;
B(2,3) = 0;
B(3,1) = - cw * c^2 * s / 6;
B(3,2) = 0;
B(3,3) = - Mbdot * c^3 * s / 8;

C(1,1) = 0;
C(1,2) = aw * c * s / 8;
C(1,3) = ac * c * s / 6;
C(2,1) = 0;
C(2,2) = - bw * c^2 * s / 6;
C(2,3) = - bc * c^2 * s / 4;
C(3,1) = 0;
C(3,2) = - cw * c^2 * s / 4;
C(3,3) = - cc * c^2 * s / 2;

% Set up for flutter solution

f = []; d = []; evalr = []; evali = []; ev = [];
icount = 0;

% Loop round for flutter solution at each velocity
% Determine first order eval solution       evals = eig(A\(rho*v^2*C + E));

for v = vstart : vinc : vend
    icount = icount + 1;
    vel(icount) = v;
    Q = [zeros(modes, modes) eye(modes, modes); -A\(rho * v^2 * C + E) -A\(rho * v * B + D)];
    evals = eig(Q);
    er = real(evals);
    ei = imag(evals);
    
    for ii = 1 : 2 * modes
        wrad(ii,1) = abs(evals(ii));
        zeta(ii,1) = -real(evals(ii)) / wrad(ii);
        whz(ii,1) = wrad(ii) / (2 * pi);
    end
    
    %   [ei, eisort] = sort(ei);    % sort evalues in order of imaginary part of evalue
    %   [er, ersort] = sort(er);    % real part of evals for plotting
    
    [whz, wsort] = sort(whz);   % sort frequencies in order
    
    f = [f whz]; d = [d zeta(wsort) * 100];
    evalr = [evalr er(wsort)]; evali = [evali ei(wsort)];
    ev = [ev evals(wsort)];
end

% Determine the flutter and divergence velocity based upon real evalue behaviour picking out crossings
% Depending upon imag evalue determine whether flutter or divergence

vs = sum(cumsum((real(ev) > 0), 2) == 1);
nvels = vs(sum(cumsum((real(ev) > 0), 2) == 1) > 0);  % number of roots that change at each stability bound
vcrit = vel(sum(cumsum((real(ev) > 0), 2) == 1) > 0);   % velocities at stability bounds

icount = 0;
for ii = 1:max(size(nvels))
    for jj = 1:nvels(ii)
        icount = icount + 1;
        v_all(icount) = vcrit(ii);     % crit velocity vector with extra (same) value for flutter root
    end
end

ecrit = ev(cumsum((real(ev) > 0), 2) == 1);

for ii = 1:max(size(ecrit))
    if imag(ecrit(ii)) == 0  % check for flutter or divergence
        [num2str(v_all(ii)) ' m/s = divergence']
    else
        fflut = abs(ecrit(ii))/(2 * pi);
        [num2str(v_all(ii)) ' m/s = flutter     '  num2str(fflut) ' Hz']
    end
end

% Plot eigenvalues

figure(2)
subplot(211)
plot(vel, evalr, '-k')
xlabel('Velocity m/s')
ylabel('Real part evals')
grid
subplot(212)
plot(vel, evali, '-k')
xlabel('Velocity m/s')
ylabel('Imag part evals')
grid

% Plot frequency and damping vs velocity
figure(3)
subplot(211)
vf = v_all(1);
xfa = [vf vf];
yf = [2 5];
plot(vel,f,'k',xfa,yf,'--k')
title(['Wing bending / torsion / control flutter - elastic axis at ',num2str(fa),'% chord'])
xlabel('Velocity (m/s)')
ylabel('Frequency (Hz)')
grid
subplot(212)
vf = v_all(1);
xfa = [vf vf];
yf = [-10 10];
y0 = [0 0];
x0 = [0 vend];
plot(vel,d,'k',xfa,yf,'--k',x0,y0,'--k')
title(['Control frequency ',  num2str(fcont,3),'Hz and hinge line at ',num2str(hline),'% chord'])
xlabel('Velocity (m/s)')
ylabel('Damping Ratio (%)')
axis([0 vend -10 15])
grid



