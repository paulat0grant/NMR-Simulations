% ====================================
% Updated twospins.m codes for 1D INEPT Spectroscopy
% Original codes for pauli.m available on 
% "https://www.oocities.org/tsm3nmr/research/nmr/mywork/mat/mat.html"
% ====================================
pauli;

IiIi=diprod(Ii,Ii);

IiIa=diprod(Ii,Ia);
IaIi=diprod(Ia,Ii);
IbIi=diprod(Ib,Ii);
IiIb=diprod(Ii,Ib);
IaIa=diprod(Ia,Ia);
IaIb=diprod(Ia,Ib);
IbIa=diprod(Ib,Ia);
IbIb=diprod(Ib,Ib);
IxIi=diprod(Ix,Ii);
IyIi=diprod(Iy,Ii);
IzIi=diprod(Iz,Ii);
IiIx=diprod(Ii,Ix);
IiIy=diprod(Ii,Iy);
IiIz=diprod(Ii,Iz);
IxIx=diprod(Ix,Ix);
IxIy=diprod(Ix,Iy);
IxIz=diprod(Ix,Iz);
IyIx=diprod(Iy,Ix);
IyIy=diprod(Iy,Iy);
IyIz=diprod(Iy,Iz);
IzIx=diprod(Iz,Ix);
IzIy=diprod(Iz,Iy);
IzIz=diprod(Iz,Iz);

IpIi=IxIi+i*IyIi;
IiIp=IiIx+i*IiIy;
ImIi=IxIi-i*IyIi;
IiIm=IiIx-i*IiIy;

Hhh=diprod(H,H);
Hih=diprod(Ii,H);
Hhi=diprod(H,Ii);

% Single transition operators:
% ----------------------------

IxIa=diprod(Ix,Ia);
IxIb=diprod(Ix,Ib);
IaIx=diprod(Ia,Ix);
IbIx=diprod(Ib,Ix);

IyIa=diprod(Iy,Ia);
IyIb=diprod(Iy,Ib);
IaIy=diprod(Ia,Iy);
IbIy=diprod(Ib,Iy);

IzIa=diprod(Iz,Ia);
IzIb=diprod(Iz,Ib);
IaIz=diprod(Ia,Iz);
IbIz=diprod(Ib,Iz);
