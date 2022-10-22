% Flutter model with control

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sets up binary aeroelastic system - flexible wing + full span control
% Continuous wing  z = q1*(y/s)^2 + q2*(y/s)*(x-xf)
% Applies PID control with Kv and Kd gains for velocity and displacement
% terms via a control surface
% Assumes that the transducer is at the wing tip leading edge.
% IAAL 2nd edition pp 208 - 212

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% IF start and end velocities are the SAME ...
%       Determines the response at one velocity from (1 - cosine) gust with and
%       without the control gains specified i.e. runs the program twice
%       Calls simulink routine   Model_H3c_ASE_Gusts_PID  (Figures 1-5)

% IF start and end velocities are DIFFERENT ...
%       Obtains flutter solution for the no gain and for the gains specified
%       at all velocities in the range specified (Figures 6-8)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format compact; clear all; close all

% Input data

% Gains  - to be used on the second loop  - first loop = zero gain

Kv =  -0.02;      % gain for velocity terms
Kd = 0;           % gain for displacement terms

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Velocity values - dictates gust or stability cases

vstart = 80; vinc = 1; vend = 140;
nvel = floor((vend - vstart) / vinc) + 1;

%%%%%%%%%   If nvel = 1 then determines the gust response for +/- gain
%%%%%%%%%   If nvel > 1 then determines stability for +/- gain

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% System parameters
s = 7.5;            % semi span
c = 2;              % chord
xcm = 0.5 * c;      % position of centre of mass from leading edge ***** DEFAULT VALUE - CAN'T BE CHANGED
xf = 0.48 * c;      % position of elastic axis from leading edge
e = xf / c - 0.25;  % eccentricity between elastic axis and aero centre (1/4 chord)
aw = 2 * pi;        % lift curve slope
bw = e * aw;        % control surface lift curve slope
m = 100;            % unit mass / area of wing
rho = 1.225;        % air density

Mtdot = -1.2;       % unsteady torsional damping term

modes = 2;

bending_freq = 5;   % bending freq in Hz - approximate - ignores coupling term
torsion_freq = 10;  % torsion freq in Hz - approximate - ignores coupling term

% Inertia matrix A
A = zeros(modes,modes);
A(1,1) = (m * s * c) / 5;
A(2,2) = m * s / 3 * (c^3 / 3 - c^2 * xf + xf^2 * c);
A(1,2) = m * s / 4 * (c^2 / 2 - c * xf);
A(2,1) = A(1,2);

% Stiffness matrix E
E = zeros(modes,modes);
EI = (bending_freq * pi * 2)^2 * A(1,1) / 4 * s^3;   % q1 term - bending stiffness
GJ = (torsion_freq * pi * 2)^2 * A(2,2) * s;         % q2 term - torsion stiffness
E(1,1) = 4 * EI / s^3;
E(2,2) = GJ / s;

EE = 0.1;     % fraction of chord made up by control surface
ac = aw / pi * (acos(1 - 2 * EE) + 2 * sqrt(EE * (1 - EE)));
bc = -aw / pi * (1 - EE) * sqrt(EE * (1 - EE));

% Aerodynamic damping matrix
B = zeros(modes,modes);
B = [c * s * aw / 10, 0; -c^2 * s * bw / 8, -c^3 * s * Mtdot / 24];

% Aerodynamic stiffness matrix
C = zeros(modes,modes);
C = [0, c * s * aw / 8; 0, -c^2 * s * bw / 6];

% (1 - cosine) gust parameters
dt = 0.001; tmin = 0; tmax = 5;
t = [0:dt:tmax]';              % Column vector
npts = max(size(t));

gust_amp_1_minus_cos = 10;     % Max velocity of (1 - cosine) gust   (m/s)
gust_t = 0.05;                 % Fraction of total time that is gust (0 - 1)

% Variables
% Variables to be updated each time round the velocity loop
freqs = []; damps = []; rpart = []; ipart = [];
% Applied control variables to be updated each time round the velocity loop
cfreqs = []; cdamps = []; crpart = []; cipart = []; vels = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Velocity loop - nvel = 1 gust response and nvel = 0 stability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for V = vstart:vinc:vend
    
    vels = [vels V];
    
    for iloop = 0:1:1 ;    %   Loop with and w/o control gain
        
        % Control surface terms (rhs of aeroelastic equations)
        F_control = rho * V^2 * [-c * s * ac / 6; c^2 * s * bc / 4];
        
        % Control feedback matrices
        g1 = F_control(1);
        g2 = F_control(2);
        F = Kv * [g1 -g1 * xf;  g2 -g2 * xf];  % control damping matrix
        G = Kd * [g1 -g1 * xf;  g2 -g2 * xf];  % control stiffness matrix
        
        % Gust terms (rhs of aeroelastic equations)
        F_gust = rho * V * [-aw * s * c / 6; bw * c^2 * s / 4];
        
        % Loop for no control and control
        if iloop == 0       % First loop without any control terms
            CC = rho * V * B;
            KK = rho * V^2 * C + E ;
        else                % Second loop is for the control term case
            CC = CC - F;
            KK = KK - G ;
        end
        
        %  Set up system matrices (in standard not aeroelastic notation)
        MM = A;
        MCC = MM\CC;
        MKK = MM\KK;
        MFG = MM\F_gust;
        
        if nvel > 1    %  Determine open / closed loop STABILITY
            
            peval = polyeig(MKK,MCC,eye(2,2));  % Solve using 2nd order form
            ff = abs(peval) / 2 / pi;
            dd = -real(peval)./abs(peval) * 100;
            rp = real(peval); ip = imag(peval);
            ff = ff(1:2:3,:); dd = dd(1:2:3,:);
            rp = rp(1:2:3,:); ip = ip(1:2:3,:);
            
            if iloop == 0       % First loop without any control terms
                freqs = [freqs ff]; damps = [damps dd];
                rpart = [rpart rp]; ipart = [ipart ip];
            else                % Second loop is for the control term case
                cfreqs = [cfreqs ff]; cdamps = [cdamps dd];
                crpart = [crpart rp]; cipart = [cipart ip];
            end
            
        else    % Otherwise compute the open / closed loop GUST RESPONSE
            
            % Gust input terms for a (1 - cosine) gust
            Sgust = zeros(size(t));
            g_end = tmax * gust_t;
            gt = sum(t < g_end);
            
            if gust_amp_1_minus_cos ~= 0
                for ii = 1:gt
                    Sgust(ii) = gust_amp_1_minus_cos / 2 * (1 - cos(2 * pi * t(ii) / g_end));
                end
            end
            
            Egust = [t,Sgust]; % Gust Array composed of time and data columns
            
            [tout] = sim('Aero_Servo');
            
            q1 = EoutG1(:,1);   % q1 generalised coordinate
            q2 = EoutG1(:,2);   % q2 generalised coordinate
            q1dot = EoutG2(:,1);
            q2dot = EoutG2(:,2);
            
            % Determine the control response (in rads)
            beta = Kv * (q1dot * s^2 - q2dot * s * xf) + Kd * (q1dot * s^2 - q2dot * s * xf);
            
            % Determine the measurement position response
            z = q1 * s^2 - q2 * s * xf;
            
            if iloop == 0
                z_nogain = z;
                q1_nogain = q1;
                q2_nogain = q2;
            end
            
        end  % evals or response loop
        
    end  % iloop loop
    
end  %  velocity loop

% Sort frequencies and corresponding dampings and poles

[freqs,iord] = sort(freqs);
nsize = size(freqs,2);
for j = 1:nsize, damps(:,j) = damps(iord(:,j),j); end
for j = 1:nsize, rpart(:,j) = rpart(iord(:,j),j); end
for j = 1:nsize, ipart(:,j) = ipart(iord(:,j),j); end

[cfreqs,iord] = sort(cfreqs);
for j = 1:nsize, cdamps(:,j) = cdamps(iord(:,j),j); end
for j = 1:nsize, crpart(:,j) = crpart(iord(:,j),j); end
for j = 1:nsize, cipart(:,j) = cipart(iord(:,j),j); end

if nvel == 1    %%% Plot out the responses for a single case
    figure(1)
    plot(t,Sgust)
    title('Gust input time history (1 - cosine)')
    xlabel('Time (s)')
    ylabel('Gust velocity (m/s)')
    
    figure(2)
    subplot(211)
    plot(t,q1,'r',t,q1_nogain,'k')
    title('q1 and q2 responses for open and closed loop cases')
    xlabel('Time (s)')
    ylabel('q1 Response ')
    grid
    legend('Closed loop','Open loop')
    subplot(212)
    plot(t,q2,'r',t,q2_nogain,'k')
    xlabel('Time (s)')
    ylabel('q2 Response ')
    grid
    
    figure(3)
    plot(t,beta)
    title(['Control Angle for Closed Loop Case','   Kv = ',num2str(Kv),...
        ' Kd = ',num2str(Kd), '  Vel = ',num2str(V), 'm/s'] )
    xlabel('Time (s)')
    ylabel('Control Angle (deg)')
    grid
    
    figure(4)
    subplot(311)
    plot(t,z_nogain,'k')
    title('Wing tip Leading Edge Displacement (m)')
    xlabel(['Kv = ',num2str(Kv),'  Kd = ',num2str(Kd), '  Vel = ',num2str(V),...
        'm/s'])
    ylabel('Open loop')
    subplot(312)
    plot(t,z,'r')
    xlabel('Time (s)')
    ylabel('Closed loop')
    subplot(313)
    plot(t,z,'r',t,z_nogain,'k')
    xlabel('Time (s)')
    ylabel('Both')
    grid
    legend('Closed loop','Open loop')
    
    figure(5)
    plot(t,z,'r',t,z_nogain,'k')
    xlabel('Time (s)')
    ylabel('Open and closed loop response')
    grid
    legend('Closed loop','Open loop')
    title(['Kv = ',num2str(Kv),'  Kd = ',num2str(Kd), '  Vel = ',num2str(V),...
        'm/s'])
    
    resp_no_gain_max = max(abs(z_nogain));
    resp_gain_max = max(abs(z));
    reduction = 100 * (resp_no_gain_max  - resp_gain_max) / resp_no_gain_max;
    astring = ['Reduction = ' num2str(reduction),  ' %'];
    disp(astring)
    
else               %%%     vg, vw, eigenvalue and root locus plots
    
    figure(6)
    subplot(211)
    plot(vels,freqs(1,:),'k',vels,freqs(2,:),'b',vels,cfreqs(1,:),'g',vels,cfreqs(2,:),'r')
    xlabel('Velocity m/s')
    ylabel('Frequency (Hz)')
    title(['Kv = ',num2str(Kv),'  Kd = ',num2str(Kd)])
    grid
    subplot(212)
    plot(vels,damps(1,:),'k',vels,damps(2,:),'b',vels,cdamps(1,:),'g',vels,cdamps(2,:),'r')
    xlabel('Velocity m/s')
    ylabel('Damping (%)')
    grid
    legend('Mode 1 open loop','Mode 2 open loop','Mode 1 closed loop','Mode 2 closed loop',...
        'Location','NorthWest')
    
    figure(7)
    subplot(211)
    plot(vels,rpart(1,:),'k',vels,rpart(2,:),'b',vels,crpart(1,:),'g',vels,crpart(2,:),'r')
    xlabel('Velocity m/s')
    ylabel('Real part of evals')
    title(['Kv = ',num2str(Kv),'  Kd = ',num2str(Kd)])
    grid
    legend('Mode 1 open loop','Mode 2 open loop','Mode 1 closed loop','Mode 2 closed loop',...
        'Location','SouthWest')
    subplot(212)
    plot(vels,ipart(1,:),'k',vels,ipart(2,:),'b',vels,cipart(1,:),'g',vels,cipart(2,:),'r')
    xlabel('Velocity m/s')
    ylabel('Imag part of evals')
    grid
    
    figure(8)
    plot(rpart(1,:),ipart(1,:),'k',rpart(2,:),ipart(2,:),'b', ...
        crpart(1,:),cipart(1,:),'g',crpart(2,:),cipart(2,:),'r',...
        rpart(1,1),ipart(1,1),'Ok',crpart(1,1),cipart(1,1),'Og',...
        rpart(2,1),ipart(2,1),'Ob',crpart(2,1),cipart(2,1),'Or')
    xlabel('Real part')
    ylabel('Imag part')
    title(['Kv = ',num2str(Kv),'  Kd = ',num2str(Kd), '  Vstart = ',num2str(vstart),...
        'm/s  Vend = ',num2str(vend), 'm/s'])
    grid
    legend('Mode 1 open loop','Mode 2 open loop','Mode 1 closed loop','Mode 2 closed loop',...
        'Location','NorthWest')
end


