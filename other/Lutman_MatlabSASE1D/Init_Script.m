alfvenCurrent = 17045.0;
mc2 = 510.99906E-3;
c   = 2.99792458E8;
e   = 1.60217733E-19;
h   = 4.13566751691*10^-15; % [eV s]
%iopt=5;
step_ratio      = T1{1};

npart			= T1{3};
energy			= T2{1};

eSpread			= T2{2};
emitN			= T2{3};
currentMax		= T2{4};
beta			= T2{5};
unduPeriod		= T3{1};
unduK			= T3{4};
% en_chirp        = IN.ebeam_bothund.energychirp;

gamma0  = energy/mc2;

if(T5{1}>=0) %Center of the spectrum (seed energy or radiation central energys)
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
% z_steps=unduL*T1{1};
coopLength = resWavelength/(4*pi*rho);			% cooperation length 
gainLength = unduPeriod/(4*pi*rho);				% rough gain length    

% z0    = unduL/gainLength;						% wiggler length in units of gain length

delt  = 4*pi*rho/unduPeriod/step_ratio; %z0/z_steps;								% integration step in z0 ~ 0.1 gain length
dels  = delt;									% integration step in s0 must be same as in z0

Ns    = currentMax/unduPeriod/step_ratio*resWavelength/c/e;
s_step = radWavelength / unduPeriod / step_ratio ; % This is the step in s watch out may be res wL and not rad...!

% gbar  = -(resWavelength-radWavelength)/(radWavelength*2*rho);		% scaled detune parameter
% delg  = eSpread/rho;							% Gaussian energy spread in units of rho
% Ns    = currentMax*unduL/unduPeriod/z_steps*resWavelength/c/e;	% N electrons per s-slice [ ]
% taper = Egain/energy/rho;        % convert Egain to taper parameter
% 
% dels2 = radWavelength * 10^6 / T1{1} / T3{1};
% 
% energy_profile=energy + ((1:s_steps)*dels*coopLength-mean((1:s_steps)*dels*coopLength))*en_chirp/bunchLength/2;
% gammas=energy_profile/mc2;
% resWavelengths=unduPeriod*(1+unduK^2/2.0)./(2*gammas.^2);
% gbars=-(resWavelengths-radWavelength)/(radWavelength*2*rho);

% s_step_size = handles.Radiation_Wavelength * 10^6 / REC.Undulator_Z_steps_ratio / REC.unduPeriod;
% (1+unduK^2/2.0)/(2*gamma0^2)/IN.general.Undulator_Z_steps_ratio;

% % % Seed Setup
% % ar = zeros(s_steps,z_steps+1);
% % ai = zeros(s_steps,z_steps+1);
% % if(strcmp(IN.starting_field.seedfile,'none'))
% %     ar(:,1) = real(input_field)'/sqrt(rhoPbeam);
% %     ai(:,1) = imag(input_field)'/sqrt(rhoPbeam);
% % else % load seed setup and overrides other paramenters (Amplitude + Phase)
% %     
% % end
    

% bunchShape		= IN.ebeam_bothund.bunchShape;
% bunchLength		= IN.ebeam_bothund.bunchLength;

% P0              = IN.starting_field.P0;
% input_field		= IN.starting_field.input_field;
% constseed		= IN.starting_field.constseed;

% IN.firstsec.z_steps=round(IN.general.Undulator_Z_steps_ratio *IN.firstsec.unduL);

% z_steps			= IN.firstsec.z_steps;
% unduL			= IN.firstsec.unduL;
% Egain			= IN.firstsec.Egain;