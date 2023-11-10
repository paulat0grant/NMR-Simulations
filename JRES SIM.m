% SIMULATION OF A THREE-SPIN 2D JRESOLVED 1H NMR SPECTRUM
% (pi/2)-(t1/2)-pi-(t1/2) sequence
% Updated codes for Jresolved Spectroscopy
% Original codes available on "https://www.oocities.org/tsm3nmr/research/nmr/mywork/mat/mat.html"
% =========================================================

clear
clc
tic
thrspins;


% INPUT PARAMETERS
% ----------------
v1=100;     % Resonance frequency of 1st spin in Hz
v2=200;     % Resonance frequency of 2nd spin in Hz
v3=300;     % Resonance frequency of 3nd spin in Hz
J12=20;     % coupling in Hz
J13=10;     % coupling in Hz
J23=15;      % coupling in Hz
T2=.6;      % Transverse decay constant in second

swh1=50;  % Spectral width for F1 Dimension
swh2=1000;  % Spectral width for F2 Dimension
td1=256;    % t1 Time domain size
td2=256;    % t2 Time domain size


% CALCULATING TIME AND FREQUENCY AXIS
% -----------------------------------
dw1=1/(swh1);                       % Dwell time
t1=0:dw1:(td1-1)*dw1;             % Time axis 
f1=(-td1/2:(td1-1)/2)/(td1*dw1);  % Frequency axis

dw2=1/swh2;                       % Dwell time
t2=0:dw2:(td2-1)*dw2;             % Time axis 
f2=(-td2/2:(td2-1)/2)/(td2*dw2);  % Frequency axis


% INITIAL DENSITY MATRIX AND HAMILTONIANS
% ---------------------------------------
indm=Izii+Iizi+Iiiz;                             % Equilibrium density matrix I1z+I2z+I3z;
csham=2*pi*(v1*Izii+v2*Iizi+v3*Iiiz);            % Chemical Shift Hamiltonian
jham=pi*J12*2*(Ixxi+Iyyi+Izzi)+pi*J13*2*(Ixix+Iyiy+Iziz)+...
     pi*J23*2*(Iixx+Iiyy+Iizz);                  % Scalar coupling Hamiltonian
eham=csham+jham;                                 % Full Hamiltonian for evolution
detop=(Ipii+Iipi+Iiip);                          % Detection Operator

% DENSITY MATRIX AFTER FIRST PULSE (pi/2)y
% ----------------------------------------
Urf=expm(-i*(Iyii+Iiyi+Iiiy)*pi/2);
outdm1=Ixii+Iixi+Iiix;

% Density matrix after 2nd Pulse pi
% ----------------------------------
Urf2 = expm(-i*(Iyii+Iiyi+Iiiy)*pi);


% EVOLUTION, MIXING, DETECTION
% ----------------------------
for m1=1:td1                        % t2 loop
   for m2=1:td2                     % t1 loop
      outdm2=expm(-i*eham*t1(m1)/2)*outdm1*expm(i*eham*t1(m1)/2)*exp(-t1(m1)/2*T2);
      outdm3=Urf2*outdm2*inv(Urf2);   % Mixing
      outdm4=expm(-i*eham*t1(m1)/2)*outdm3*expm(i*eham*t1(m1)/2)*exp(-t1(m1)/2*T2);
      outdm5=expm(-i*eham*t2(m2))*outdm4*expm(i*eham*t2(m2))*exp(-t2(m2)/T2);
      s(m1,m2)=trace(detop*outdm5); % Detection
   end
end    


% FOURIER TRANSFORM AND PLOTTING
% ------------------------------
S=fftshift(fft2(-s)); % Fourier Transform

contour(f2,f1,real(S),20);
toc
    
