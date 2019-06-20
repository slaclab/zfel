z_steps = round(UndulatorLength*T1{1});
taper = Egain/energy/rho;

% figure, plot(abs(Final_field)), title('final (at start)')

ar = zeros(s_steps,z_steps+1);
ai = zeros(s_steps,z_steps+1);

ar(:,1)=real(Final_field(:,1));
ai(:,1)=imag(Final_field(:,1));

% figure, plot(ar(:,1))

thethalf = zeros(npart,z_steps+1);
gam = zeros(npart,z_steps+1);

final_theta=CurrentPhaseSpace.theta;
final_gamma=CurrentPhaseSpace.gamma;

for k = 1:s_steps-1
        if shape(k)>0.05                                       % calculate FEL interaction
            thet0=CurrentPhaseSpace.theta(k,:);
            gam0=CurrentPhaseSpace.gamma(k,:);
%             [thet0,gam0] = load_bucket(npart,initial_gamma(k),deviation_gamma(k),iopt,Ns);	% load each bucket
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
%                 if(Save_phase_space_n_steps>0)
%                     if(mod(j-1,Save_phase_space_n_steps)==0)
%                         dump_particle_inserted=dump_particle_inserted+1;
%                         if(~IN.saving.Save_phase_space_only_avg_std )
%                         dump_particle_gamma_distance(k,:,dump_particle_inserted)=gam(:,j+1);
%                         dump_particle_theta_distance(k,:,dump_particle_inserted)=thethalf(:,j+1);
%                         else
%                            avg_gamma(k,dump_particle_inserted)=mean(gam(:,j+1));
%                            std_gamma(k,dump_particle_inserted)=std(gam(:,j+1));
%                         end
%                     end                    
%                 end
            end
            final_gamma(k,:)=gam(:,z_steps+1);
            final_theta(k,:)=thethalf(:,z_steps+1);
        else
            final_gamma(k,:)=gam0; % ignore FEL interaction
            final_theta(k,:)=thet0;
            for j = 1:z_steps                     
                ar(k+1,j+1) = ar(k,j);              % apply slippage condition
                ai(k+1,j+1) = ai(k,j);      
            end
        end
end

Final_field=ar(:,end)+1i*ai(:,end);

% figure, plot(abs(Final_field)), title('final')

CurrentPhaseSpace.theta=final_theta;
CurrentPhaseSpace.gamma=final_gamma;