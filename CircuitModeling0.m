% PA 7 for ELEC 4700

close all
clear 
clc

% v = i*r;
% i = c*dv/dt;
% V = ldi*dt;

% Circuit Components
r1 = 1;
g1 = 1/r1;
r2 = 2;
g2 = 1/r2;
r3 = 22; % Obtained with Assigment 3 code
g3 = 1/r3;
r4 = 0.1;
g4 = 1/r4;
ro = 1000;
go = 1/ro;

c = 0.25;
L = 0.2;
a = 100;
vin = 1;

% eq1: derivative(v1-v2)/dt*C + (v1-v2)*g1 + Is = 0;
% eq2: v1 = vin;
% eq3: derivative(v2-v1)/dt*C + (v2-v1)*g1 + v2*g2 + Il = 0;
% eq4: v3*g3 - Il = 0;
% eq5: i3 = v3*v4;
% eq6: (v4-v5)*g4 + i3 = 0;
% eq7: (v5-v4)*g4 + v5*g0 = 0; 
% eq8: v3-v2 = L*dIl/dt;

% creating matrices
X = ['v1';'v2';'Is';'v3';'i3';'v4';'v5';'Il'];

G = [g1,-g1,1,0,0,0,0,0;
    1,0,0,0,0,0,0,0;
    -g1,g1+g2,0,0,0,0,0,1;
    0,0,0,g3,0,0,0,-1;
    0,0,0,0,1,0,0,0;
    0,0,0,0,1,g4,-g4,0;
    0,0,0,0,0,-g4,g4+go,0;
    0,-1,0,1,0,0,0,0];
    
C = [c,-c,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0;
    -c,c,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,L];

% F = [0,vin,0,0,0,0,0];
F = zeros(1,length(X));
F(1,2) = vin;

% CdV/dt + GV = F
% (G + jwC)X = F(w)
% F(w) = 1
w = pi;
H = G+1j*w*C;
X = H\F.';

% DC sweep
w = 0;
vin = linspace(-10,10,100);

for v = vin
    H = G+1j*w*C;
    F = zeros(1,length(X));
    F(1,2) = v;
    
    X0 = H\F.';
    v3 = 5*X0(8,:);
    vo = 5*X0(4,:);

    figure(1)
    plot(v,vo,'r.')
    hold on
    grid on
    plot(v,v3,'b.')
    legend('V_o','V_3')
    title('DC Voltage Plots')
    xlabel('V_{in} (V)'),ylabel('Voltage (V)')
end

% AC sweep
vin = 1;
F(1,2) = vin;
w = 1:1:100;
F(1,2) = vin;

for w0 = w
    
    H = G+w0*C;
    X0 = H\F.';
    v3 = 5*X0(8,:);
    vo = 10*X0(4,:);
    
    figure(2) 
    plot(w0,vo,'r.')
    hold on
    grid on
    plot(w0,v3,'b.')
    legend('V_o','V_3')
    title('AC Voltage Plots')
    xlabel('{\omega} (Hz)'),ylabel('Voltage (V)')
end

% Capacitance Sweeps
n = 1000;
mu = 0.25;
sigma = 0.05;
c0 = zeros(1,n);
gain = zeros(1,n);
w0 = pi;
vin = 1;

for k = 1:n
    c0(k) = normrnd(mu,sigma);
    C = [c0(k),-c0(k),0,0,0,0,0,0;
    0,0,0,0,0,0,0,0;
    -c0(k),c0(k),0,0,0,0,0,0;
    0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,L];

    figure(3)
    histogram(c0,20)
    title('Capacitance Distribution')
    xlabel('C'),ylabel('Count')
    
    H = G+w0*c0(k);
    X1 = F/H;
    vo(k) = X1(:,7);
    gain(k) = vo(k)/vin;

    figure(4) 
    histogram(gain,20)
    title('Voltage Gain Distribution')
    xlabel('V_o/V_{in}'),ylabel('Count')    
end

