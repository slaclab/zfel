alfvenCurrent = 17045.0;
mc2 = 510.99906E-3;
c   = 2.99792458E8;
e   = 1.60217733E-19;
h   = 4.13566751691*10^-15; % [eV s]
step_ratio      = T1{1};
npart			= T1{3};
energy			= T2{1};
eSpread			= T2{2};
emitN			= T2{3};
currentMax		= T2{4};
beta			= T2{5};
unduPeriod		= T3{1};
unduK			= T3{4};
gamma0  = energy/mc2;

if(T5{1}>0) %Center of the spectrum (seed energy or radiation central energys)
    radWavelength	= T5{2}*10^-9;
else
    radWavelength   = unduPeriod*(1+unduK^2/2)/2/gamma0^2;
end

if(T1{4})
    iopt=5;
else
    iopt=4;
end

s_steps			= T1{2};
unduJJ  = besselj(0,unduK^2/(4+2*unduK^2))-besselj(1,unduK^2/(4+2*unduK^2));
gamma0  = energy/mc2;							% mean gamma
sigmaX2 = emitN*beta/gamma0;					% sigmaX square
rho     = (0.5/gamma0)*((currentMax/alfvenCurrent)*(unduPeriod*unduK*unduJJ/(2*pi))^2/(2*sigmaX2))^(1/3);
resWavelength = unduPeriod*(1+unduK^2/2.0)/(2*gamma0^2);  % resonant wavelength to gamma0 THIS IS FOR ELECTRON DETUNING !!!
rhoPbeam   = rho*energy*1e6*currentMax;
coopLength = resWavelength/(4*pi*rho);			% cooperation length 
gainLength = unduPeriod/(4*pi*rho);				% rough gain length    
delt  = 4*pi*rho/unduPeriod/step_ratio; %z0/z_steps;								% integration step in z0 ~ 0.1 gain length
dels  = delt;									% integration step in s0 must be same as in z0
Ns    = currentMax/unduPeriod/step_ratio*resWavelength/c/e;
s_step = radWavelength / unduPeriod / step_ratio ; % This is the step in s watch out may be res wL and not rad...!
s_step_safe=s_step;
step_in_z=1/T1{1};