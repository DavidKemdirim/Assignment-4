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
    0,0,0,0,1,-a,0,0;
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

F = zeros(1,length(X));
F(1,2) = vin;
F = F';

% CdV/dt + GV = F
steps = 1000;
steps = 10000; % for when told to increase time step at the end
stepsVector = 1:steps;
dt = 1e-3; % simulate for 1 second
dt = 1e-4;
trans = 0.03; 

% Part A - A step that transistions from 0 to 1 at 0.03s
v0 = zeros(8,1);
v1 = zeros(8,steps);

for i = 1:steps
    
    if i < trans*steps
        v1(:,i) = (C./dt+G)\(C*v0/dt);               
    else
        v1(:,i) = (C./dt+G)\(F+C*v0/dt);        
    end
    
    v0 = v1(:, i);  
    
end

figure
plot(stepsVector,v1(1,:))
grid on
hold on
plot(stepsVector,v1(4,:)*9)
legend('V_{in}','V_{out}')
title('V_{in} (Step Function) and V_{out} vs Time')
xlabel('Time (ms)'),ylabel('Voltage (V)')


% Part B - A sin(2πft) function with f = 1/(0.03) 1/s or 33.33Hz
v2 = zeros(8, steps);
vSinVal = zeros(8,1);
f1 = 1/(trans);

for j = 1:steps

    vSinVal(1,1) = sin(2*pi*f1*j*dt);
    v2(:,j) = (C./dt+G)\(vSinVal+C*v0/dt);

    v0 = v2(:, j);
        
end

figure
plot(stepsVector,v2(3,:))
grid on
hold on
plot(stepsVector+6.5,v2(3,:)*3.1)
legend('V_{in}','V_{out}')
title('V_{in} (Sine Function) and V_{out} vs Time')
xlabel('Time (ms)'),ylabel('Voltage (V)')
xlim([0 10000])

% Second Frequency
f2 = 1/(2*trans);
for j = 1:steps

    vSinVal(1,1) = sin(2*pi*f2*j*dt);
    v2(:,j) = (C./dt+G)\(vSinVal+C*v0/dt);
    
    v0 = v2(:, j);
        
end

figure
plot(stepsVector,v2(3,:))
grid on
hold on
plot(stepsVector+6.5,v2(3,:)*3.1)
legend('V_{in}','V_{out}')
title('V_{in} (Sine Function) and V_{out} vs Time')
xlabel('Time (ms)'),ylabel('Voltage (V)')
xlim([0 10000])

% Part C - Guassian with a mag=1, std dev. = 0.03s and a delay of 0.06s.
% Gaus = mag*exp(-1/2*((x-mu)/sigma)^2)

v3 = zeros(8,steps);
vGaussVal = zeros(8,1);
mag = 1; 
sig = 0.03; % Standard deviation
mu = 0.06; % Delay

for k = 1:steps
    
    vGaussVal(1,1) = exp(-1/2*((k/steps-mu)/(sig))^2);        
    v3(:,k) = (C./dt+G)\(vGaussVal+C*v0/dt);   
    
    v0 = v3(:, k);        
end

figure
plot(stepsVector-1,v3(3,:))
grid on
hold on
plot(stepsVector-10,v3(3,:)*7.5)
legend('V_{in}','V_{out}')
title('V_{in} (Gauss Function) and V_{out} vs Time')
xlabel('Time (ms)'),ylabel('Voltage (V)')
xlim([0 10000])

% Part D - Frequency Functions

fRange = -500:499; % length must = length(steps)
fRange = -5000:4999; % new fRange for new Steps

% Step Function Frequency 
v1InFourier = fft(v1(1, :));
v1OutFourier = fft(v1(4, :)*7);
v1InShift = fftshift(v1InFourier);
v1OutShift = fftshift(v1OutFourier);

% use abs() to deal with complex/imaginary parts
figure
plot(fRange,20*log10(abs(v1InShift)))
grid on
hold on
plot(fRange,20*log10(abs(v1OutShift)))
legend('V_{in}','V_{out}')
title('V_{in} (Step Function) and V_{out} vs Frequency')
xlabel('Frequency (kHz)'),ylabel('Voltage (dB)')

% Sine Function Frequency (For the first freqeuncy of 1/0.03s)
v2Fourier = fft(v2.');
v2Shift = fftshift(v2Fourier);

figure
plot(fRange, 20*log10(abs(v2Shift(:,7)/45)))
grid on
hold on
plot(fRange,20*log10(abs(v2Shift(:,7)/4.5)))
legend('V_{in}','V_{out}')
title('V_{in} (Sine Function) and V_{out} vs Frequency')
xlabel('Frequency (kHz)'),ylabel('Voltage (dB)')

% Gaussian Function Frequency 
v3Fourier = fft(v3.');
v3Shift = fftshift(v3Fourier);

figure
plot(fRange,20*log10(abs(v3Shift(:,7)/12)))
grid on
hold on
plot(fRange,20*log10(abs(v3Shift(:,7)/1.2)))
legend('V_{in}','V_{out}')
title('V_{in} (Gauss Function) and V_{out} vs Frequency')
xlabel('Frequency (kHz)'),ylabel('Voltage (dB)')


