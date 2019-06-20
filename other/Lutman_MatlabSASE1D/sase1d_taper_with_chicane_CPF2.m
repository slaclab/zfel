function [OUT]=sase1d_taper_with_chicane_CPF2(IN)

npart			= IN.ebeam_bothund.npart;
energy			= IN.ebeam_bothund.energy;

eSpread			= IN.ebeam_bothund.eSpread;
emitN			= IN.ebeam_bothund.emitN;
currentMax		= IN.ebeam_bothund.currentMax;
beta			= IN.ebeam_bothund.beta;
unduPeriod		= IN.ebeam_bothund.unduPeriod;
unduK			= IN.ebeam_bothund.unduK;
en_chirp        = IN.ebeam_bothund.energychirp;
fixed_s_axis    = IN.ebeam_bothund.fix_energy_for_s_at;

radWavelength	= IN.ebeam_bothund.radWavelength;
s_steps			= IN.general.s_steps;
bunchShape		= IN.ebeam_bothund.bunchShape;
bunchLength		= IN.ebeam_bothund.bunchLength;

P0              = IN.starting_field.P0;
input_field		= IN.starting_field.input_field;
constseed		= IN.starting_field.constseed;

z_steps			= IN.firstsec.z_steps;
unduL			= IN.firstsec.unduL;
Egain			= IN.firstsec.Egain;

Save_phase_space_n_steps=IN.saving.Save_phase_space_n_steps;
Save_field_every_n_steps=IN.saving.Save_field_every_n_steps;
BothFieldandPower=IN.saving.BothFieldandPower;
SeedInformation=IN.saving.SeedInformation;

Use_Average_gamma_and_deviation=IN.chicaneandcrystal.Use_Average_gamma_and_deviation;
T0=IN.chicaneandcrystal.T0;                      
WakeLength=IN.chicaneandcrystal.WakeLength;                   
ElectronDelay=IN.chicaneandcrystal.ElectronDelay;               

alfvenCurrent = 17045.0;
mc2 = 510.99906E-3;
c   = 2.99792458E8;
e   = 1.60217733E-19;
h   = 4.13566751691*10^-15; % [eV s]
iopt=5;

unduJJ  = besselj(0,unduK^2/(4+2*unduK^2))-besselj(1,unduK^2/(4+2*unduK^2));
gamma0  = energy/mc2;							% mean gamma
sigmaX2 = emitN*beta/gamma0;					% sigmaX square
rho     = (0.5/gamma0)*((currentMax/alfvenCurrent)*(unduPeriod*unduK*unduJJ/(2*pi))^2/(2*sigmaX2))^(1/3);
resWavelength = unduPeriod*(1+unduK^2/2.0)/(2*gamma0^2);  % resonant wavelength to gamma0
rhoPbeam   = rho*energy*1e6*currentMax;         % rho times beam power [W]
coopLength = resWavelength/(4*pi*rho);			% cooperation length 
gainLength = unduPeriod/(4*pi*rho);				% rough gain length    
z0    = unduL/gainLength;						% wiggler length in units of gain length
delt  = z0/z_steps;								% integration step in z0 ~ 0.1 gain length
dels  = delt;									% integration step in s0 must be same as in z0
gbar  = -(resWavelength-radWavelength)/(radWavelength*2*rho);		% scaled detune parameter
delg  = eSpread/rho;							% Gaussian energy spread in units of rho
Ns    = currentMax*unduL/unduPeriod/z_steps*resWavelength/c/e;	% N electrons per s-slice [ ]
taper = Egain/energy/rho;        % convert Egain to taper parameter

energy_profile=energy + ((1:s_steps)*dels*coopLength-mean((1:s_steps)*dels*coopLength))*en_chirp/bunchLength/2;
gammas=energy_profile/mc2;
resWavelengths=unduPeriod*(1+unduK^2/2.0)./(2*gammas.^2);
gbars=-(resWavelengths-radWavelength)/(radWavelength*2*rho);

% Electron Bunch Setup
if(strcmp(IN.ebeam_bothund.ebeam_bunchfile,'none'))
    bunch_steps=round(bunchLength/delt/coopLength);  %rms (Gaussian) or half width (flattop) bunch length in s_step
    shape = zeros(1,s_steps);
    switch(bunchShape)
        case 1
            shape= 0.5*(tanh(10*((1:s_steps)-s_steps/2+bunch_steps)/bunch_steps)-tanh(10*((1:s_steps)-s_steps/2-bunch_steps)/bunch_steps));
%             load GutsoX CurrentProfile
%             shape=CurrentProfile/3000;
        case 2
            shape= exp(-((1:s_steps)-s_steps/2).^2/(2*bunch_steps^2));
        case 3
            shape=shape+1;
        case 4
            shape(1:(bunch_steps*2))=1;
        case 5
            shape= 0.5*(tanh(10*((1:s_steps)-s_steps/2+bunch_steps)/bunch_steps)-tanh(10*((1:s_steps)-s_steps/2-bunch_steps)/bunch_steps))/4; 
            shape=shape+exp(-((1:6000)-4500).^2/10000)*0.75;
            figure, plot(shape)
        case 6
            shape= 0.5*(tanh(10*((1:s_steps)-s_steps/2+bunch_steps)/bunch_steps)-tanh(10*((1:s_steps)-s_steps/2-bunch_steps)/bunch_steps));
            load GutsoX CurrentProfile
            shape=CurrentProfile/3000;
    end
    
    shape(10000:11800)=shape(10000:11800)*0.9;
    
    initial_gamma=gbars;
    deviation_gamma=ones(size(shape))*delg;
else % load bunch setup and overrides other paramenters (Shape + energy + energy deviation)
    eval(['load ',IN.ebeam_bothund.ebeam_bunchfile])
    currentMax=max(CurrentProfile);
    energy=CurrentProfile*EnergyProfile.'/(sum(CurrentProfile));
    gamma0  = energy/mc2;
    delt  = z0/z_steps;								
    dels  = delt;
    sstepsize = (1+unduK^2/2.0)/(2*gamma0^2)/IN.general.Undulator_Z_steps_ratio;
    sigmaX2 = emitN*beta/gamma0;
    rho     = (0.5/gamma0)*((currentMax/alfvenCurrent)*(unduPeriod*unduK*unduJJ/(2*pi))^2/(2*sigmaX2))^(1/3);
    resWavelength = unduPeriod*(1+unduK^2/2.0)/(2*gamma0^2);
    rhoPbeam   = rho*energy*1e6*currentMax;
    coopLength = resWavelength/(4*pi*rho);
    gainLength = unduPeriod/(4*pi*rho);
    z0    = unduL/gainLength;
    gbar  = -(resWavelength-radWavelength)/(radWavelength*2*rho);
    delg  = eSpread/rho;
    Ns    = currentMax*unduL/unduPeriod/z_steps*resWavelength/c/e;
    taper = Egain/energy/rho;
    s_steps=round((max(sProfile)-min(sProfile))/sstepsize/10^6)+IN.firstsec.z_steps+IN.secondsec.z_steps;
    s_steps_no_space=round((max(sProfile)-min(sProfile))/sstepsize/10^6);
    
    energy_profile=interp1(sProfile,EnergyProfile,linspace(min(sProfile),max(sProfile),s_steps_no_space));
    current_profile=interp1(sProfile,CurrentProfile,linspace(min(sProfile),max(sProfile),s_steps_no_space));
    espread_profile=interp1(sProfile,espreadProfile,linspace(min(sProfile),max(sProfile),s_steps_no_space));
    
    energy_profile=[energy_profile,zeros(1,IN.firstsec.z_steps+IN.secondsec.z_steps)];
    current_profile=[current_profile,zeros(1,IN.firstsec.z_steps+IN.secondsec.z_steps)];
    espread_profile=[espread_profile,zeros(1,IN.firstsec.z_steps+IN.secondsec.z_steps)];
    
    gammas=energy_profile/mc2;
    resWavelengths=unduPeriod*(1+unduK^2/2.0)./(2*gammas.^2);
    gbars=-(resWavelengths-radWavelength)/(radWavelength*2*rho);
    
    initial_gamma=gbars;
    deviation_gamma=espread_profile/rho;
    shape=current_profile/max(current_profile);
    input_field=zeros(size(shape));
end

% Seed Setup
ar = zeros(s_steps,z_steps+1);
ai = zeros(s_steps,z_steps+1);
if(strcmp(IN.starting_field.seedfile,'none'))
    ar(:,1) = real(input_field)'/sqrt(rhoPbeam);
    ai(:,1) = imag(input_field)'/sqrt(rhoPbeam);
else % load seed setup and overrides other paramenters (Amplitude + Phase)
    
end

% Array spaces setup
final_gamma = zeros(s_steps,npart);
final_theta = zeros(s_steps,npart);
if(Save_phase_space_n_steps>0)
    if(~IN.saving.Save_phase_space_only_avg_std )
        dump_particle_gamma_distance = zeros(s_steps,npart,sum(mod(0:z_steps,Save_phase_space_n_steps)==0));
        dump_particle_theta_distance = zeros(s_steps,npart,sum(mod(0:z_steps,Save_phase_space_n_steps)==0));
    else
        avg_gamma=zeros(s_steps,sum(mod(0:z_steps,Save_phase_space_n_steps)==0));
        std_gamma=zeros(s_steps,sum(mod(0:z_steps,Save_phase_space_n_steps)==0));
    end
end
thethalf = zeros(npart,z_steps+1);
gam = zeros(npart,z_steps+1);

%random seed control
if constseed==1
        rand('state',0);						% hold constant random seed if requested
        randn('state',0);                       % hold constant random seed if requested
end

% FIRST SECTION  
    
 % go over all slices of the bunch starting from the tail k=1 for shape>0.1
    for k = 1:s_steps-1
        if shape(k)>0.05                                       % calculate FEL interaction
            [thet0,gam0] = load_bucket(npart,initial_gamma(k),deviation_gamma(k),iopt,Ns);	% load each bucket
            gam(:,1) = gam0';							% gamma at j=1
            thethalf(:,1) = thet0'-gam(:,1)*delt/2;     % half back
            dump_particle_inserted=0;
            for j = 1:z_steps                           % evolve e and gamma in s and t by leap-frog
                thet = thethalf(:,j)+gam(:,j)*delt/2+taper(j)*delt/2;
                sumsin = sum(sin(thet));
                sumcos = sum(cos(thet));
                sinavg = shape(k)*sumsin/npart;
                cosavg = shape(k)*sumcos/npart;
                arhalf = ar(k,j)+cosavg*dels/2;
                aihalf = ai(k,j)-sinavg*dels/2;
                thethalf(:,j+1) = thethalf(:,j)+gam(:,j)*delt+taper(j)*delt;
                gam(:,j+1) = gam(:,j)-2*arhalf*cos(thethalf(:,j+1))*delt...
                    +2*aihalf*sin(thethalf(:,j+1))*delt;
                sumsin = sum(sin(thethalf(:,j+1)));
                sumcos = sum(cos(thethalf(:,j+1)));
                sinavg = shape(k)*sumsin/npart;
                cosavg = shape(k)*sumcos/npart;
                ar(k+1,j+1) = ar(k,j)+cosavg*dels;	% apply slippage condition
                ai(k+1,j+1) = ai(k,j)-sinavg*dels;
                if(Save_phase_space_n_steps>0)
                    if(mod(j-1,Save_phase_space_n_steps)==0)
                        dump_particle_inserted=dump_particle_inserted+1;
                        if(~IN.saving.Save_phase_space_only_avg_std )
                        dump_particle_gamma_distance(k,:,dump_particle_inserted)=gam(:,j+1);
                        dump_particle_theta_distance(k,:,dump_particle_inserted)=thethalf(:,j+1);
                        else
                           avg_gamma(k,dump_particle_inserted)=mean(gam(:,j+1));
                           std_gamma(k,dump_particle_inserted)=std(gam(:,j+1));
                        end
                    end                    
                end
            end
            final_gamma(k,:)=gam(:,z_steps+1);
            final_theta(k,:)=thethalf(:,z_steps+1);
        else
            final_gamma(k,:)=zeros(1,npart); % ignore FEL interaction
            final_theta(k,:)=zeros(1,npart);
            for j = 1:z_steps                     
                ar(k+1,j+1) = ar(k,j);              % apply slippage condition
                ai(k+1,j+1) = ai(k,j);      
            end
        end
    end
    
    % calculate radiated energy as a function of z and output power 
    OUT.firstsection.z=(0:1:z_steps)*delt*gainLength;
   
    if(Save_phase_space_n_steps>0)
        OUT.firstsection.dump_particle_z=(0:Save_phase_space_n_steps:z_steps)*delt*gainLength;
        if(~IN.saving.Save_phase_space_only_avg_std )
            OUT.firstsection.dump_particle_gamma_distance=dump_particle_gamma_distance;
            OUT.firstsection.dump_particle_theta_distance=dump_particle_theta_distance;
        else
            OUT.firstsection.gamma_avg=avg_gamma;
            OUT.firstsection.gamma_std=std_gamma;
        end
    end
    
    
    if(Save_phase_space_n_steps==0)
        OUT.firstsection.dump_particle_gamma_distance=final_gamma;
        OUT.firstsection.dump_particle_theta_distance=final_theta;
        OUT.firstsection.dump_particle_z=(z_steps)*delt*gainLength;
    end
    a2 = ar.^2+ai.^2;
    OUT.firstsection.power_z = mean(a2)*rhoPbeam;
    OUT.firstsection.energy_z = sum(a2)*dels*coopLength/c*rhoPbeam*1e3;
    OUT.firstsection.s = (1:s_steps)*dels*coopLength*1.0e6;                        % longitudinal coordinate [micron]
    if(~isnan(fixed_s_axis))
        OUT.firstsection.s = (1:s_steps)*1.0e6*(1+unduK^2/2.0)/(2*(fixed_s_axis/mc2)^2)/IN.general.Undulator_Z_steps_ratio;
    else
        OUT.firstsection.s = (1:s_steps)*dels*coopLength*1.0e6; 
    end
    
    if(Save_field_every_n_steps)
        OUT.firstsection.z_for_field=(0:Save_field_every_n_steps:z_steps)*delt*gainLength;
        OUT.firstsection.field_s=complex(ar(:,1:Save_field_every_n_steps:(z_steps+1)),ai(:,1:Save_field_every_n_steps:(z_steps+1)))*sqrt(rhoPbeam);
        if(BothFieldandPower)
            OUT.firstsection.power_s=abs(OUT.firstsection.field_s).^2;
        end
    end
    OUT.firstsection.detune = 2*pi/(dels*s_steps)*(-s_steps/2:s_steps/2-1);      % detune in units of 2*rho relative to radWavelength
    OUT.firstsection.exit_field = complex(ar(:,z_steps+1),ai(:,z_steps+1))*sqrt(rhoPbeam);
    if(BothFieldandPower)
            OUT.firstsection.exit_power=abs(OUT.firstsection.exit_field).^2;
    end
    OUT.firstsection.detune0=gbar;                                               % output seed frequency detune
    OUT.firstsection.rho=rho;
    OUT.firstsection.shape=shape;
    OUT.firstsection.resWavelength=resWavelength;
    OUT.firstsection.PhotonEnergyeV=h/resWavelength*c;
    OUT.firstsection.EnergyProfile=energy_profile;
    OUT.firstsection.Initial_gammas=gammas;
    clear a2;clear ar; clear ai; clear final_theta;
if(IN.general.simulate_only_firstsection ~=1)
    
 % CHICANE
 StepLength=OUT.firstsection.s(2)-OUT.firstsection.s(1);
 s1=0:StepLength:WakeLength; 
 TransmittedBraggG=(besselj(1,sqrt(s1/T0))./sqrt(s1/T0))/(2*T0);
 TransmittedBraggG(1)=1/4/T0;
 seed=conv(OUT.firstsection.exit_field(end:-1:1),TransmittedBraggG)*StepLength; 
 seed=seed(end:-1:1);
 New_s_coordinate=[s1(2:end),s1(end)+OUT.firstsection.s ];
 First_seed_position=s1(end)-ElectronDelay;
 [First_seed_position_value,First_seed_position_index]=min(abs(New_s_coordinate- First_seed_position));
 New_seed=seed(First_seed_position_index:(-1+First_seed_position_index+length(OUT.firstsection.s )));
 New_SeedToSlip=seed;
 New_SeedToSlip(First_seed_position_index:(-1+First_seed_position_index+length(OUT.firstsection.s )))=0;
 
 if(SeedInformation)
     OUT.chicane.TransmittedBraggG=TransmittedBraggG;
     OUT.chicane.s_TransmittedBraggG=s1;
     OUT.chicane.seed_full_wake=seed;
     OUT.chicane.seed_full_wake_with_bunchhole=New_SeedToSlip;
     OUT.chicane.s_seed_full_wake=New_s_coordinate;
     OUT.chicane.seed_over_the_electrons=New_seed;
     OUT.chicane.s=OUT.firstsection.s;
 end

 % Setting up for radiator part
 if(~IN.general.simulate_only_firstsection)
z_steps			= IN.secondsec.z_steps;
unduL			= IN.secondsec.unduL;
Egain			= IN.secondsec.Egain;

unduJJ  = besselj(0,unduK^2/(4+2*unduK^2))-besselj(1,unduK^2/(4+2*unduK^2));
gamma0  = energy/mc2;							% mean gamma
sigmaX2 = emitN*beta/gamma0;					% sigmaX square
rho     = (0.5/gamma0)*((currentMax/alfvenCurrent)*(unduPeriod*unduK*unduJJ/(2*pi))^2/(2*sigmaX2))^(1/3);
resWavelength = unduPeriod*(1+unduK^2/2.0)/(2*gamma0^2);  % resonant wavelength to gamma0
rhoPbeam   = rho*energy*1e6*currentMax;         % rho times beam power [W]
coopLength = resWavelength/(4*pi*rho);			% cooperation length 
gainLength = unduPeriod/(4*pi*rho);				% rough gain length
%cs0  = bunchLength/coopLength					% bunch length in units of cooperation length     
z0    = unduL/gainLength;						% wiggler length in units of gain length
delt  = z0/z_steps;								% integration step in z0 ~ 0.1 gain length
dels  = delt;									% integration step in s0 must be same as in z0
gbar  = -(resWavelength-radWavelength)/(radWavelength*2*rho);		% scaled detune parameter
delg  = eSpread/rho;							% Gaussian energy spread in units of rho
Ns    = currentMax*unduL/unduPeriod/z_steps*resWavelength/c/e;	% N electrons per s-slice [ ]
taper = Egain/energy/rho;        % convert Egain to taper parameter

%Electron bunch set up
if(Use_Average_gamma_and_deviation)
    initial_gamma=ones(size(shape))*(shape*mean(final_gamma,2))/sum(shape);
    deviation_gamma=ones(size(shape))*(shape*std(final_gamma.').')/sum(shape);
else
    initial_gamma=mean(final_gamma,2).';
    deviation_gamma=std(final_gamma.');
end

%Seed Set up
ar = zeros(s_steps,z_steps+1);
ai = zeros(s_steps,z_steps+1);
ar(:,1) = real(New_seed)'/sqrt(rhoPbeam);
ai(:,1) = imag(New_seed)'/sqrt(rhoPbeam);

% Array spaces setup
final_gamma = zeros(s_steps,npart);
final_theta = zeros(s_steps,npart);

if(Save_phase_space_n_steps>0)
    if(~IN.saving.Save_phase_space_only_avg_std )
        dump_particle_gamma_distance = zeros(s_steps,npart,sum(mod(0:z_steps,Save_phase_space_n_steps)==0));
        dump_particle_theta_distance = zeros(s_steps,npart,sum(mod(0:z_steps,Save_phase_space_n_steps)==0));
    else
        avg_gamma=zeros(s_steps,sum(mod(0:z_steps,Save_phase_space_n_steps)==0));
        std_gamma=zeros(s_steps,sum(mod(0:z_steps,Save_phase_space_n_steps)==0));
    end
end

thethalf = zeros(npart,z_steps+1);
gam = zeros(npart,z_steps+1);

%random seed control
if constseed==1
        rand('state',0);						% hold constant random seed if requested
        randn('state',0);                       % hold constant random seed if requested
end

    for k = 1:s_steps-1
        if shape(k)>0.05                                       % calculate FEL interaction
            [thet0,gam0] = load_bucket(npart,initial_gamma(k),deviation_gamma(k),iopt,Ns);	% load each bucket
            gam(:,1) = gam0';							% gamma at j=1
            thethalf(:,1) = thet0'-gam(:,1)*delt/2;     % half back
            dump_particle_inserted=0;
            for j = 1:z_steps                           % evolve e and gamma in s and t by leap-frog
                thet = thethalf(:,j)+gam(:,j)*delt/2+taper(j)*delt/2;
                sumsin = sum(sin(thet));
                sumcos = sum(cos(thet));
                sinavg = shape(k)*sumsin/npart;
                cosavg = shape(k)*sumcos/npart;
                arhalf = ar(k,j)+cosavg*dels/2;
                aihalf = ai(k,j)-sinavg*dels/2;
                thethalf(:,j+1) = thethalf(:,j)+gam(:,j)*delt+taper(j)*delt;
                gam(:,j+1) = gam(:,j)-2*arhalf*cos(thethalf(:,j+1))*delt...
                    +2*aihalf*sin(thethalf(:,j+1))*delt;
                sumsin = sum(sin(thethalf(:,j+1)));
                sumcos = sum(cos(thethalf(:,j+1)));
                sinavg = shape(k)*sumsin/npart;
                cosavg = shape(k)*sumcos/npart;
                ar(k+1,j+1) = ar(k,j)+cosavg*dels;	% apply slippage condition
                ai(k+1,j+1) = ai(k,j)-sinavg*dels;
                if(Save_phase_space_n_steps>0)
                    if(mod(j-1,Save_phase_space_n_steps)==0)
                        dump_particle_inserted=dump_particle_inserted+1;
                        if(~IN.saving.Save_phase_space_only_avg_std )
                        dump_particle_gamma_distance(k,:,dump_particle_inserted)=gam(:,j+1);
                        dump_particle_theta_distance(k,:,dump_particle_inserted)=thethalf(:,j+1);
                        else
                           avg_gamma(k,dump_particle_inserted)=mean(gam(:,j+1));
                           std_gamma(k,dump_particle_inserted)=std(gam(:,j+1));
                        end
                    end                   
                end
            end
            final_gamma(k,:)=gam(:,z_steps+1);
            final_theta(k,:)=thethalf(:,z_steps+1);
        else
            final_gamma(k,:)=zeros(1,npart); % ignore FEL interaction
            final_theta(k,:)=zeros(1,npart);
            for j = 1:z_steps                     
                ar(k+1,j+1) = ar(k,j);              % apply slippage condition
                ai(k+1,j+1) = ai(k,j);      
            end
        end
    end
    
    % calculate radiated energy as a function of z and output power 
    OUT.secondsection.z=(0:1:z_steps)*delt*gainLength;
    if(Save_phase_space_n_steps>0)
        OUT.secondsection.dump_particle_z=(0:Save_phase_space_n_steps:z_steps)*delt*gainLength;
        if(~IN.saving.Save_phase_space_only_avg_std )
            OUT.secondsection.dump_particle_gamma_distance=dump_particle_gamma_distance;
            OUT.secondsection.dump_particle_theta_distance=dump_particle_theta_distance;
        else
            OUT.secondsection.gamma_avg=avg_gamma;
            OUT.secondsection.gamma_std=std_gamma;
        end
    end
    if(Save_phase_space_n_steps==0)
        OUT.secondsection.dump_particle_gamma_distance=final_gamma;
        OUT.secondsection.dump_particle_theta_distance=final_theta;
        OUT.secondsection.dump_particle_z=(z_steps)*delt*gainLength;
    end
    a2 = ar.^2+ai.^2;
    OUT.secondsection.power_z = mean(a2)*rhoPbeam;
    OUT.secondsection.energy_z = sum(a2)*dels*coopLength/c*rhoPbeam*1e3;
    OUT.secondsection.s = (1:s_steps)*dels*coopLength*1.0e6;                        % longitudinal coordinate [micron]
    
    if(~isnan(fixed_s_axis))
        OUT.firstsection.s = (1:s_steps)*1.0e6*(1+unduK^2/2.0)/(2*(fixed_s_axis/mc2)^2)/IN.general.Undulator_Z_steps_ratio;
    else
        OUT.firstsection.s = (1:s_steps)*dels*coopLength*1.0e6; 
    end
    
    if(Save_field_every_n_steps)
        OUT.secondsection.z_for_field=(0:Save_field_every_n_steps:z_steps)*delt*gainLength;
        OUT.secondsection.field_s=complex(ar(:,1:Save_field_every_n_steps:(z_steps+1)),ai(:,1:Save_field_every_n_steps:(z_steps+1)))*sqrt(rhoPbeam);
        if(BothFieldandPower)
            OUT.secondsection.power_s=abs(OUT.secondsection.field_s).^2;
        end
    end
    OUT.secondsection.detune = 2*pi/(dels*s_steps)*(-s_steps/2:s_steps/2-1);      % detune in units of 2*rho relative to radWavelength
    OUT.secondsection.exit_field = complex(ar(:,z_steps+1),ai(:,z_steps+1))*sqrt(rhoPbeam);
    OUT.secondsection.detune0=gbar;  
    OUT.secondsection.rho=rho;
    OUT.secondsection.shape=shape;
    OUT.secondsection.resWavelength=resWavelength;
    OUT.secondsection.PhotonEnergyeV=h/resWavelength*c;
    
 end

end

if(IN.processing.CalculateSpectrum)
    ZeroPadLength=[1,1]*IN.processing.Number_of_zero_padding;
    BaseRange=OUT.firstsection.detune*2*OUT.firstsection.rho; %Linspace of the undulator period.
    Orel=linspace(BaseRange(1),BaseRange(length(BaseRange)),ZeroPadLength(1)+ZeroPadLength(2)+length(OUT.firstsection.detune));
    STDom=length(OUT.firstsection.detune);
    CNorm=(OUT.firstsection.s(2)-OUT.firstsection.s(1));
    if(Save_field_every_n_steps)
    zmax=length(OUT.firstsection.z_for_field);
    fieldFFT_z=zeros(length(Orel),zmax);
    for Distance_Counter=1:zmax
        fieldFFT_z(:,Distance_Counter) = CNorm*circshift(fft([zeros(ZeroPadLength(1),1);OUT.firstsection.field_s(:,Distance_Counter);zeros(ZeroPadLength(2),1)]),round((STDom+ZeroPadLength(1)+ZeroPadLength(2))/2));
    end
    OUT.firstsection.fieldFFT_z=fieldFFT_z;
    OUT.firstsection.Orel=Orel;
    if(BothFieldandPower)
        OUT.firstsection.Pspec_z=abs(fieldFFT_z).^2;
    end
    else
        OUT.firstsection.exit_fieldFFT_z= CNorm*circshift(fft([zeros(ZeroPadLength(1),1);OUT.firstsection.exit_field;zeros(ZeroPadLength(2),1)]),round((STDom+ZeroPadLength(1)+ZeroPadLength(2))/2));
        OUT.firstsection.Orel=Orel;
        if(BothFieldandPower)
            OUT.firstsection.exit_Pspec_z=abs(OUT.firstsection.exit_fieldFFT_z).^2;
        end
    end
    
    if(~IN.general.simulate_only_firstsection)
        BaseRange=OUT.secondsection.detune*2*OUT.secondsection.rho;
        Orel=linspace(BaseRange(1),BaseRange(length(BaseRange)),ZeroPadLength(1)+ZeroPadLength(2)+length(OUT.secondsection.detune));
        CNorm=(OUT.secondsection.s(2)-OUT.secondsection.s(1));
        STDom=length(OUT.secondsection.detune);
        if(Save_field_every_n_steps)
        zmax=length(OUT.secondsection.z_for_field);
        fieldFFT_z=zeros(length(Orel),zmax);   
        for Distance_Counter=1:zmax
            fieldFFT_z(:,Distance_Counter) = CNorm*circshift(fft([zeros(ZeroPadLength(1),1);OUT.secondsection.field_s(:,Distance_Counter);zeros(ZeroPadLength(2),1)]),round((STDom+ZeroPadLength(1)+ZeroPadLength(2))/2));
        end
        OUT.secondsection.fieldFFT_z=fieldFFT_z;
        OUT.secondsection.Orel=Orel;
        if(BothFieldandPower)
            OUT.secondsection.Pspec_z=abs(fieldFFT_z).^2;
        end
        else
            OUT.secondsection.exit_fieldFFT_z= CNorm*circshift(fft([zeros(ZeroPadLength(1),1);OUT.secondsection.exit_field;zeros(ZeroPadLength(2),1)]),round((STDom+ZeroPadLength(1)+ZeroPadLength(2))/2));
            OUT.secondsection.Orel=Orel;
            if(BothFieldandPower)
                OUT.secondsection.exit_Pspec_z=abs(OUT.secondsection.exit_fieldFFT_z).^2;
            end
        end
    end    
end
