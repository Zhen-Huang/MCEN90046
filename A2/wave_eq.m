% wave_eq.m: solve wave equation response to an initial condition
clear all, close all
l = 1; c = 1; % choose length, wave speed
N = 200; % choose number of modes for modal solution
sig = 0.03; % choose \sigma for Gaussian initial condition

% this bit turns initial condition into Gamma_i(0) using fancy methods to
% integrate -- you can ignore it
[x,quads] = clencurt(N); x = (1-x)*l/2; quads = quads*l/2;
% q0 = sin(2*pi*1*x); % one sinusoid
% q0 = sin(2*pi*5*x); % more sinusoids
% q0 = sin(2*pi*1*x) + sin(2*pi*2*x); % two sinusoids
 q0 = 0.01*exp(-(x-0.5*l).^2./sig^2); % gaussian

Q = sin(pi*x/l*[1:N]); gamma0 = zeros(N,1);
for i=1:N
  integrand = q0 .* Q(:,i); gamma0(i) = dot(quads,integrand);
end

% now to simulate using the Gamma_i values as the initial condition
t = [0:0.005:10]; x = [0:0.001:1]'*l; % choose time and grid
% find N natural frequencies and mode shapes
om = pi*l/c*[1:N]; Q = sin(pi*x/l*[1:N]);
% add in some damping to each mode; or set to zero
% zeta = 0.001 * [1:N]; 
zeta = zeros(size(om));

% set up state-space matrices
A = [zeros(N) eye(N); -diag(om.^2) -diag(2*zeta.*om)];
B = [zeros(N); eye(N)]; C = [Q zeros(length(x),N)];
H = ss(A,B,C,0); g = zeros(length(t),N); % forcing is zero
q = lsim(H,g,t,[gamma0;zeros(size(gamma0))]); q=q';
%
% for i = 1:3:length(t)
%   fig=figure(1); clf; hold on
%   set(gca,'Box','On','XMinorTick','On','YMinorTick','On');         
%   plot(x,q(:,1),'bo-'), plot(x,q(:,i),'ko-'), axis([0 1 -.6 .6]),
%   xlabel('t'), ylabel('q(x,t)'), title(sprintf('t=%f',t(i)))
% end

q = q*1000;
figure(2)
hold on;
plot(x,q(:,1));
plot(x,q(:,4));
plot(x,q(:,8));
plot(x,q(:,20));
plot(x,q(:,40));
legend('0 us','20 us','40 us','100 us','200 us')
xlabel('L/m')
ylabel('q(x,t)/mm')
