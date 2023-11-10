% ====================================
% SIMULATION OF A TWO-SPIN 1D SPECTRUM
% T. S. Mahesh (20 May 2002)
% ====================================

% NOTE: This program uses diprod.m, pauli.m and twospins.m.
% All these files should be in the same folder.

clear
clc

twospins

% INPUT PARAMETERS
% ----------------
v1=-100;    % Resonance frequency of 1st spin in Hz
v2=200;     % Resonance frequency of 2nd spin in Hz
J=30;       % coupling in Hz
T2=.5;      % Transverse decay constant in second

swh=1000;   % Spectral width in Hz
td=4*1024;  % Time domain size
theta1=pi/2; % Pulse angle
theta2=pi    % Pulse angle during pi
JIS = 145;  % Coupling constant between I and S

% CALCULATING TIME AND FREQUENCY AXIS
% -----------------------------------
dw=1/swh;           % Dwell time
t=0:dw:(td-1)*dw;   % Time axis 
k=-td/2:(td-1)/2;
f=k/(td*dw);        % Frequency axis
d2=1/(4*JIS)        % Delay during spin echo


% INITIAL DENSITY MATRIX AND HAMILTONIANS
% ---------------------------------------
indm=IzIi+IiIz;                  % Equilibrium density matrix I1z+I2z;
csham=2*pi*(v1*IzIi+v2*IiIz);    % Chemical Shift Hamiltonian
jham=pi*J*2*(IxIx+IyIy+IzIz);     % Scalar coupling Hamiltonian
eham=csham+jham;               % Full Hamiltonian for evolution
detop=(IpIi+IiIp);               % Detection Operator

% DENSITY MATRIX AFTER RF
% -----------------------
Urf1=expm(-i*(IyIi+IiIy)*theta1);  % For theta(y) pulse
Urf2=expm(-i*(IyIi+IiIy)*theta2);   % For theta(y) pulse during echo
outdm1=Urf1*indm*inv(Urf1);      % Output density matrix after

% EVOLUTION AND FID
% -----------------
for m=1:td
    outdm2=expm(-i*eham*d2)*outdm1*expm(i*eham*d2);      % first delay period
    outdm3=Urf2*outdm2*inv(Urf2);                         % 180 pulse
    outdm4=expm(-i*eham*d2)*outdm3*expm(i*eham*d2);      % 2nd delay period
    outdm5=Urf1*outdm4*inv(Urf1);                        % 90 pulse
    outdm6=expm(-i*eham*d2)*outdm5*expm(i*eham*d2);      % 3rd delay period
    outdm7=Urf2*outdm6*inv(Urf2);                         % 180 pulse
    outdm8=expm(-i*eham*d2)*outdm7*expm(i*eham*d2);      % 4th delay period
    outdm9=expm(-i*eham*t(m))*outdm8*expm(i*eham*t(m));  % t2 evolution
    s(m)=trace(detop*outdm9);                            % Detection
end    
s=s.*exp(-t/T2); 					 % T2 decay

% FOURIER TRANSFORM AND PLOTTING
% ------------------------------
S=fftshift(fft(s)); % Fourier Transform
subplot(2,1,1),plot(t,real(s));
xlabel('Time in second'); ylabel('Intensity (s)');
title('Simulation of a two-spin system');
subplot(2,1,2),plot(f,real(S));
xlabel('Frequency in Hz'); ylabel('Intensity (S)');
