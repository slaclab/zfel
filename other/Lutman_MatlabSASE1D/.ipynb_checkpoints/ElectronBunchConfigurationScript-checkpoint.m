%bunch_steps=round(bunchLength/delt/coopLength);  %rms (Gaussian) or half width (flattop) bunch length in s_step
shape = zeros(1,T1{2});
pulse_duration=T2{6}*c/10^15;

TotalWidth=s_step*T1{2};
    
    switch(BunchShape)
        case 1 %flat top hard edges
            bunch_steps=round(pulse_duration/s_step_safe);
            shape(1:min(bunch_steps,T1{2}))=1;
        case 2 %flat top soft edges
            bunch_steps=round(pulse_duration/s_step_safe)/2;
            shape= 0.5*(tanh(10*((1:s_steps)-s_steps/2+bunch_steps)/bunch_steps)-tanh(10*((1:s_steps)-s_steps/2-bunch_steps)/bunch_steps));
        case 3 %Gaussian
            sigmagaussian=pulse_duration/sqrt(log(256))/s_step;
            shape=exp(-((1:T1{2}) - T1{2}/2 ).^2/2/(sigmagaussian^2));
            shape=shape/max(shape);
        case 4 %from file
            load(NameShape)
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

gamma0  = energy/mc2;    
resWavelength = unduPeriod*(1+unduK^2/2.0)/(2*gamma0^2);
gbar  = -(resWavelength-radWavelength)/(radWavelength*2*rho);		% scaled detune parameter
delg  = eSpread/rho;                                                % Gaussian energy spread in units of rho
%taper = Egain/energy/rho;        % convert Egain to taper parameter
initial_gamma=ones(size(shape))*gbar;
deviation_gamma=ones(size(shape))*delg;

CurrentPhaseSpace.theta=zeros(length(shape),npart);
CurrentPhaseSpace.gamma=zeros(length(shape),npart);

for k = 1:s_steps-1
    [thet0,gam0] = load_bucket(npart,initial_gamma(k),deviation_gamma(k),5-T1{4},Ns);	% load each bucket
    CurrentPhaseSpace.theta(k,:)=thet0;
    CurrentPhaseSpace.gamma(k,:)=gam0;
%     gam(:,1) = gam0';							% gamma at j=1
%     thethalf(:,1) = thet0'-gam(:,1)*delt/2; 
end
