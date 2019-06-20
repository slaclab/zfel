% General Simulation Parameters sampling, beam parameters ... INPUT AREA
IN.general.Undulator_Z_steps_ratio  = 2;        % integration steps = Undulator_Z_steps_ratio*length_meters +1
IN.general.Number_of_Shots          = 1;        % Number of shots to simulate per Run
IN.general.s_steps = 6000;                    % samples on longitudinal for the simulation
IN.general.simulate_only_firstsection = 1;      % 1 simulation stops before chicane ; 2 also simulates chicane ; 0 simulates also second section
IN.general.Number_of_Runs          = 1;        % Number of Runs
IN.general.Save_On_Disk          = 0;           % 1 saves on disk 0 no saves.

IN.ebeam_bothund.npart               = 1024; % n-macro-particles per simulated slice
IN.ebeam_bothund.meanenergy          = 5.810595384217249e+03; % electron energy [MeV] center
IN.ebeam_bothund.energy_deviation_rms = 0;  % each shot has a random energy using this rms from center
IN.ebeam_bothund.eSpread = 0.5e-4;  %uncorrelated slice energy spread
IN.ebeam_bothund.emitN   = 0.4e-6;		% normalized transverse emittance [mm-mrad]
IN.ebeam_bothund.currentMax = 3000.0;					% peak current [Ampere]
IN.ebeam_bothund.beta = 25.0;							% mean beta [meter]
IN.ebeam_bothund.unduPeriod = 0.03;						% undulator period [meter]
IN.ebeam_bothund.unduK = 3.5;							% undulator parameter, K [ ]
IN.ebeam_bothund.radWavelength =  8.265612869101367e-10;     		% seed wavelength [meter] or SASE spectral plot center wavelength 
IN.ebeam_bothund.bunchShape =1;                         % 1 is smoothed flattop, 2 is Gaussian, 3 is Constant, 4 real flat top
IN.ebeam_bothund.bunchLength = 2.05e-6*1.5;                    % half width for flattop, rms for gaussian [m] 
IN.ebeam_bothund.ebeam_bunchfile='none';                % a filename overrides other information
IN.ebeam_bothund.scan_energy=0;                         % [FROM,TO] scan energy linearly, 0 does not scan energy
IN.ebeam_bothund.energychirp=0;                         % Head-Tail for the bunch length
IN.ebeam_bothund.fix_energy_for_s_at=5.810595384217249e+03;

IN.processing.CalculateSpectrum = 1;           % 1 calculates spectrum when calculates the field; 0 doesn't.
IN.processing.Number_of_zero_padding = 0;

IN.saving.Save_field_every_n_steps = 1; % save full beam data every n steps. 0 Saves only the end of the simulation
IN.saving.Save_phase_space_n_steps = 1; % save full beam data every n steps. 0 Saves only the end of the simulation -1 Does not save phase space
IN.saving.Save_phase_space_only_avg_std = 1; %Do not save full particle phase space but only average and energy deviation
IN.saving.BothFieldandPower =        0;           % 1 retains both field and power, 0 just the field
IN.saving.SeedInformation =          1;  % 1 saves  information related to chicane, like seed, wake, ...

IN.firstsec.unduL = 80;							% length of undulator [meter]
IN.firstsec.zstart=200;                             % starting point of taper
IN.firstsec.alpha= 280.0e-3;                        % taper coefficient

IN.secondsec.unduL = 1;							% length of undulator [meter]
IN.secondsec.zstart=200;                             % starting point of taper
IN.secondsec.alpha= 280.0e-3;                        % taper coefficient

IN.starting_field.P0=0;                 % small seed input power [W]
IN.starting_field.input_field = sqrt(IN.starting_field.P0)*exp(-((1:IN.general.s_steps)-IN.general.s_steps/2).^2/(2*(5*IN.general.s_steps/2)^2));  % scaled initial input field vector [sqrt(W)] 
IN.starting_field.constseed=0;
IN.starting_field.seedfile='none';              % a filename overrides other information

IN.chicaneandcrystal.Use_Average_gamma_and_deviation=0; % Uses the average value for gamma and deviation.
IN.chicaneandcrystal.WakeLength=200;                         %WakeLength in Micron
IN.chicaneandcrystal.ElectronDelay=6.2;                     %Delay in Micron
IN.chicaneandcrystal.Reflection=[4,0,0];
IN.chicaneandcrystal.thickness=104;                         % crystal thickness in micron
IN.chicaneandcrystal.Bragg='b';                              % 'b' = Bragg ; 'l' =Laue
IN.chicaneandcrystal.Eta=NaN;                                % NaN= symmetric (either Bragg or Laue), otherwise will be used as Asymmetry angle.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AREA
%h   = 4.13566751691*10^-15; % [eV s]
%c   = 299792458; [m/s]
%Energy = h/resWavelength*c [eV]

[IN.chicaneandcrystal.T0,IN.chicaneandcrystal.Angle]=CrystalT0(IN.chicaneandcrystal.Reflection,4.13566751691*10^-15/(IN.ebeam_bothund.radWavelength)*299792458 ,IN.chicaneandcrystal.thickness,IN.chicaneandcrystal.Bragg,IN.chicaneandcrystal.Eta);
IN.chicaneandcrystal.T0=IN.chicaneandcrystal.T0*299792458*1e6; %in micron
IN.firstsec.z_steps=round(IN.general.Undulator_Z_steps_ratio *IN.firstsec.unduL);
IN.secondsec.z_steps=round(IN.general.Undulator_Z_steps_ratio *IN.secondsec.unduL);

IN.firstsec.Egain= ((1:IN.firstsec.z_steps) >= IN.firstsec.zstart).*((1:IN.firstsec.z_steps)*IN.firstsec.unduL - IN.firstsec.zstart).^2*IN.firstsec.alpha;
IN.secondsec.Egain=((1:IN.secondsec.z_steps) >= IN.secondsec.zstart).*((1:IN.secondsec.z_steps)*IN.secondsec.unduL - IN.secondsec.zstart).^2*IN.secondsec.alpha;

IN.secondsec.alpha=0;
IN.secondsec.zstart=22*IN.general.Undulator_Z_steps_ratio;
IN.firstsec.Egain= ((1:IN.firstsec.z_steps) >= IN.firstsec.zstart).*((1:IN.firstsec.z_steps) - IN.firstsec.zstart).^2*IN.firstsec.alpha;
IN.secondsec.Egain=((1:IN.secondsec.z_steps) >= IN.secondsec.zstart).*((1:IN.secondsec.z_steps)- IN.secondsec.zstart).^2*IN.secondsec.alpha;

for JL=1:IN.general.Number_of_Runs

for IK=1:IN.general.Number_of_Shots

    disp(['Simulating Shot ',num2str(IK+(JL-1)*IN.general.Number_of_Shots)])
    if(IN.ebeam_bothund.scan_energy)
        IN.ebeam_bothund.energy=IN.ebeam_bothund.scan_energy(1) + (IK-1)*diff(IN.ebeam_bothund.scan_energy)/(IN.general.Number_of_Shots-1);
    else
        IN.ebeam_bothund.energy=IN.ebeam_bothund.meanenergy * (1+ randn(1)*IN.ebeam_bothund.energy_deviation_rms);
    end
    OUTPUT{IK}=sase1d_taper_with_chicane_CPF(IN);
    INPUT{IK}=IN;
end
    clear IK
    if(IN.general.Save_On_Disk)
        eval(['save TempoRex0_',num2str(JL),' -v7.3'])
        clear INPUT
        clear OUTPUT
    end
end

%Make Output for the profile reconstruction code
KEEPDISTANCE=22;

OUTPUT{1}.firstsection

