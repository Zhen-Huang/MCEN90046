clear;
m = 200;
s = 7.5;
c = 2;
xf = 0.48*c;
aw = 2*pi;
e = (xf - 0.25 * c)/c;
M_th = -1.2;
rho = 1.225;
EI = 2*10^7;
GJ = 2*10^6;

A = m.*[s*c/5 s*c/6 (s/4)*((c^2)/2-c*xf);
    s*c/6 s*c/7 (s/5)*((c^2)/2-c*xf)
    (s/4)*((c^2)/2-c*xf) (s/5)*((c^2)/2-c*xf) (s/3)*((c^3)/3-c^2*xf+c*xf^2)];
B = [c*s*aw/10 c*s*aw/12 0;
    c*s*aw/12 c*s*aw/14 0;
    -c^2*s*e*aw/8 -c^2*s*e*aw/10 -c^3*s*M_th/24];
C = [0  0 c*s*aw/8;
    0  0 c*s*aw/10;
    0 0 -c^2*s*e*aw/6];
D = zeros(3,3);
E = [4*EI/(s^3) 0 0;
    0 4*EI/(s^3) 0
    0 0 GJ/s];
%% 10.3
A = m.*[s*c/5 s*c/6 (s/4)*((c^2)/2-c*xf);
    s*c/6 s*c/7 (s/5)*((c^2)/2-c*xf)
    (s/4)*((c^2)/2-c*xf) (s/5)*((c^2)/2-c*xf) (s/3)*((c^3)/3-c^2*xf+c*xf^2)];
B = zeros(3,3);
C = [0  0 c*s*aw/8;
    0  0 c*s*aw/10;
    0 0 -c^2*s*e*aw/6];
D = zeros(3,3);
E = [4*EI/(s^3) 0 0;
    0 4*EI/(s^3) 0
    0 0 GJ/s];
f =[];
d =[];
Vstart = 0;
Vinc = 0.1;
Vend = 180;
for V = Vstart:Vinc:Vend
    Q = [zeros(3,3) eye(3);
        -A\(rho*V^2*C+E) -A\(rho*V*B+D)];
    eigvalue = eig(Q);
    for ii = 1:6
        wrad(ii,1) = abs(eigvalue(ii));
        zeta(ii,1) = -real(eigvalue(ii)) / wrad(ii);
        if imag(eigvalue(ii)) == 0
            wrad(ii,1) = 0;
        end
        whz(ii,1) = wrad(ii,1) / (2*pi);
    end
    [whz, wsort] = sort(whz);
    zeta = zeta(wsort);
    
    f = [f whz]; d = [d zeta*100];
end
vel = Vstart:Vinc:Vend;
figure(1)
subplot(2,1,1);
plot(vel,f)
title("Figure 10.3 Q2.2")
xlabel('Velocity m/s'); ylabel('Frequecy hz'); grid on
subplot(2,1,2);
plot(vel,d,'.','MarkerSize',1);
xlabel('Velocity m/s'); ylabel('Damping ratio %');grid on

%% 10.4
M_th = 0;
A = m.*[s*c/5 s*c/6 (s/4)*((c^2)/2-c*xf);
    s*c/6 s*c/7 (s/5)*((c^2)/2-c*xf)
    (s/4)*((c^2)/2-c*xf) (s/5)*((c^2)/2-c*xf) (s/3)*((c^3)/3-c^2*xf+c*xf^2)];
B = [c*s*aw/10 c*s*aw/12 0;
    c*s*aw/12 c*s*aw/14 0;
    -c^2*s*e*aw/8 -c^2*s*e*aw/10 -c^3*s*M_th/24];
C = [0  0 c*s*aw/8;
    0  0 c*s*aw/10;
    0 0 -c^2*s*e*aw/6];
D = zeros(3,3);
E = [4*EI/(s^3) 0 0;
    0 4*EI/(s^3) 0
    0 0 GJ/s];
f =[];
d =[];
Vstart = 10;
Vinc = 0.1;
Vend = 90;
for V = Vstart:Vinc:Vend
    Q = [zeros(3,3) eye(3);
        -A\(rho*V^2*C+E) -A\(rho*V*B+D)];
    eigvalue = eig(Q);
    for ii = 1:6
        wrad(ii,1) = abs(eigvalue(ii));
        zeta(ii,1) = -real(eigvalue(ii)) / wrad(ii);
        if imag(eigvalue(ii)) == 0
            wrad(ii,1) = 0;
        end
        whz(ii,1) = wrad(ii,1) / (2*pi);
    end
    [whz, wsort] = sort(whz);
    zeta = zeta(wsort);
    
    f = [f whz]; d = [d zeta*100];
end
vel = Vstart:Vinc:Vend;
figure(2)
subplot(2,1,1);
plot(vel,f)
title("Figure 10.4 Q2.2")
xlabel('Velocity m/s'); ylabel('Frequecy hz'); grid on
subplot(2,1,2);
plot(vel,d,'.','MarkerSize',1);
xlabel('Velocity m/s'); ylabel('Damping ratio %');grid on
%% 10.5
clear;
m = 200;
s = 7.5;
c = 2;
aw = 2*pi;
rho = 1.225;
EI = 2*10^7;
GJ = 2*10^6;

vflu = [];
r = 0.25:0.05:0.65;
M_th = 0;
for k = 1:length(r)
    xf = r(k)*c;
    e = (xf - 0.25 * c)/c;
    A = m.*[s*c/5 s*c/6 (s/4)*((c^2)/2-c*xf);
        s*c/6 s*c/7 (s/5)*((c^2)/2-c*xf)
        (s/4)*((c^2)/2-c*xf) (s/5)*((c^2)/2-c*xf) (s/3)*((c^3)/3-c^2*xf+c*xf^2)];
    B = [c*s*aw/10 c*s*aw/12 0;
        c*s*aw/12 c*s*aw/14 0;
        -c^2*s*e*aw/8 -c^2*s*e*aw/10 -c^3*s*M_th/24];
    C = [0  0 c*s*aw/8;
        0  0 c*s*aw/10;
        0 0 -c^2*s*e*aw/6];
    D = zeros(3,3);
    E = [4*EI/(s^3) 0 0;
        0 4*EI/(s^3) 0
        0 0 GJ/s];
    f =[];
    d =[];
    ev = [];
    Vstart = 10;
    Vinc = 0.1;
    Vend = 300;
    for V = Vstart:Vinc:Vend
        Q = [zeros(3,3) eye(3);
            -A\(rho*V^2*C+E) -A\(rho*V*B+D)];
        eigvalue = eig(Q);
        for ii = 1:6
            wrad(ii,1) = abs(eigvalue(ii));
            zeta(ii,1) = -real(eigvalue(ii)) / wrad(ii);
            if imag(eigvalue(ii)) == 0
                wrad(ii,1) = 0;
            end
            whz(ii,1) = wrad(ii,1) / (2*pi);
        end
        [whz, wsort] = sort(whz);
        zeta = zeta(wsort);
        
        f = [f whz]; d = [d zeta*100];
        ev = [ev eigvalue(wsort)];
    end
    vel = Vstart:Vinc:Vend;
    vs = sum(cumsum((real(ev) > 0),2) == 1);
    nvels = vs(sum(cumsum((real(ev) > 0),2) == 1) > 0); 	% number of roots that change at each stability bound
    vcrit = vel(sum(cumsum((real(ev) > 0),2) == 1) > 0);    % velocities at stability bounds
    
    if isempty(vcrit)
        vflu = [vflu 0];
    else
        vflu = [vflu vcrit(1)];
    end
end

figure(3)
grid on
plot(r,vflu,'-^');
title("Figure 10.5 Q2.2")
hold on

vflu = [];
r = 0.25:0.05:0.65;
M_th = -1.2;
for k = 1:length(r)
    xf = r(k)*c;
    e = (xf - 0.25 * c)/c;
    A = m.*[s*c/5 s*c/6 (s/4)*((c^2)/2-c*xf);
        s*c/6 s*c/7 (s/5)*((c^2)/2-c*xf)
        (s/4)*((c^2)/2-c*xf) (s/5)*((c^2)/2-c*xf) (s/3)*((c^3)/3-c^2*xf+c*xf^2)];
    B = [c*s*aw/10 c*s*aw/12 0;
        c*s*aw/12 c*s*aw/14 0;
        -c^2*s*e*aw/8 -c^2*s*e*aw/10 -c^3*s*M_th/24];
    C = [0  0 c*s*aw/8;
        0  0 c*s*aw/10;
        0 0 -c^2*s*e*aw/6];
    D = zeros(3,3);
    E = [4*EI/(s^3) 0 0;
        0 4*EI/(s^3) 0
        0 0 GJ/s];
    f =[];
    d =[];
    ev = [];
    Vstart = 1;
    Vinc = 0.1;
    Vend = 300;
    for V = Vstart:Vinc:Vend
        Q = [zeros(3,3) eye(3);
            -A\(rho*V^2*C+E) -A\(rho*V*B+D)];
        eigvalue = eig(Q);
        for ii = 1:6
            wrad(ii,1) = abs(eigvalue(ii));
            zeta(ii,1) = -real(eigvalue(ii)) / wrad(ii);
            if imag(eigvalue(ii)) == 0
                wrad(ii,1) = 0;
            end
            whz(ii,1) = wrad(ii,1) / (2*pi);
        end
        [whz, wsort] = sort(whz);
        zeta = zeta(wsort);
        f = [f whz]; d = [d zeta*100];
        ev = [ev eigvalue(wsort)];
    end
    vel = Vstart:Vinc:Vend;
    vs = sum(cumsum((real(ev) > 0),2) == 1);
    nvels = vs(sum(cumsum((real(ev) > 0),2) == 1) > 0);	% number of roots that change at each stability bound
    vcrit = vel(sum(cumsum((real(ev) > 0),2) == 1) > 0);    % velocities at stability bounds
    
    if isempty(vcrit)
        vflu = [vflu 0];
    else
        vflu = [vflu vcrit(1)];
    end
end
plot(r,vflu,'-s');
legend("M_ theta_dot=0","M_ theta_dot=-1.2")
xlabel("Elastic axis position/chord");ylabel("Flutter speed (m/s)");

%% 10.6
clear;
m = 200;
s = 7.5;
c = 2;
xf = 0.48*c;
aw = 2*pi;
e = (xf - 0.25 * c)/c;
rho = 1.225;
EI = 2*10^7;
GJ = 2*10^6;
M_th = -1.2;

A = m.*[s*c/5 s*c/6 (s/4)*((c^2)/2-c*xf);
     s*c/6 s*c/7 (s/5)*((c^2)/2-c*xf)
    (s/4)*((c^2)/2-c*xf) (s/5)*((c^2)/2-c*xf) (s/3)*((c^3)/3-c^2*xf+c*xf^2)];
B = [c*s*aw/10 c*s*aw/12 0;
    c*s*aw/12 c*s*aw/14 0;
    -c^2*s*e*aw/8 -c^2*s*e*aw/10 -c^3*s*M_th/24];
C = [0  0 c*s*aw/8;
    0  0 c*s*aw/10;
    0 0 -c^2*s*e*aw/6];
D = zeros(3,3);
E = [4*EI/(s^3) 0 0;
    0 4*EI/(s^3) 0
    0 0 GJ/s];
f =[];
d =[];
Vstart = 10;
Vinc = 0.1;
Vend = 100;
for V = Vstart:Vinc:Vend
    Q = [zeros(3,3) eye(3);
        -A\(rho*V^2*C+E) -A\(rho*V*B+D)];
    eigvalue = eig(Q);
    for ii = 1:6
        wrad(ii,1) = abs(eigvalue(ii));
        zeta(ii,1) = -real(eigvalue(ii)) / wrad(ii);
        if imag(eigvalue(ii)) == 0
            wrad(ii,1) = 0;
        end
        whz(ii,1) = wrad(ii,1) / (2*pi);
    end
    [whz, wsort] = sort(whz);
    zeta = zeta(wsort);
    
    f = [f whz]; d = [d zeta*100];
end
vel = Vstart:Vinc:Vend;
figure(4)
subplot(2,1,1);
plot(vel,f(1:2,:),'-^','MarkerIndices',[1 101 201 301 401 501 601 701 801 901]);
title("Figure 10.6 Q2.2")
hold on
plot(vel,f(3:4,:),'-s','MarkerIndices',[1 101 201 301 401 501 601 701 801 901]);
plot(vel,f(5:6,:),'-o','MarkerIndices',[1 101 201 301 401 501 601 701 801 901]);
xlabel('Velocity m/s'); ylabel('Frequecy hz'); grid on
subplot(2,1,2);
plot(vel,d(1:2,:),'-^','MarkerIndices',[1 101 201 301 401 501 601 701 801 901]);
hold on
plot(vel,d(3:4,:),'-s','MarkerIndices',[1 101 201 301 401 501 601 701 801 901]);
plot(vel,d(5:6,:),'-o','MarkerIndices',[1 101 201 301 401 501 601 701 801 901]);
xlabel('Velocity m/s'); ylabel('Damping ratio %');grid on