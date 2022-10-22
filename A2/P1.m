bl = 0:0.02:30;
plot(bl,cos(bl).*cosh(bl)-1)
%%
bl = [4 8 11 14 17 20];
bl = fsolve(@(bl) cos(bl).*cosh(bl)-1, bl, optimset('Display','off'))
om = bl.^2;
n = length(om);

l = 1;
x = [0:0.001:1] * l;
Q = zeros(length(x),n);
for i = 1:n
    B = (sin(bl(i))-sinh(bl(i)))/(cosh(bl(i))-cos(bl(i)));
    Q(:,i) = sin(bl(i)*x)+sinh(bl(i)*x)+B*(cos(bl(i)*x)+cosh(bl(i)*x));
end

figure
for i = 1:n
    subplot(3,2,i),plot(x,Q(:,i),'k');
end

%%
L2h = (22.37/(2*pi*32.25))*sqrt(200*10^9/(12*7900));
L = sqrt(L2h*0.015);
w = om(1:6).*sqrt(200*10^9/(12*7900))/L2h;
f = w/(2*pi)