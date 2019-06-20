function [thet,gam] = load_bucket(n,gbar,delg,iopt,Ns)

%c       Load thet & gam (longitudinal phase space)
%c       iopt = 4  = > uniform in theta, Fawley beamlet for seeded
%c       iopt = 5  = > Shot noise start (Penman algorithm, Fawley beamlet)

nmax = 1000000;
if(n>nmax)
  error('increase nmax, subr load')
end
if(iopt==4)
  M=4;  % number of particles in each beamlet
  nb= round(n/M); % number of beamlet via Fawley between 64 to 256 (x16=1024 to 4096)
  if(M*nb~=n)
	error('n must be a multiple of 4')
  end

% make hammersley pseudo-random sequence of size 1xnb
  dim_num = 1; step = 0; seed = 1; leap = 1; base = 3;
  ham_rand = i_to_hammersley_sequence(dim_num, nb, step, seed, leap, base);
 
  for i=1:nb
    gamma=delg*randn(1)+gbar;
      for j=1:M
         gam((i-1)*M+j)=gamma;
         thet((i-1)*M+j)=2*pi*j/M + 2*pi*ham_rand(i);
      end
  end

% older loading without Hammersley
%   for i=1:nb
%     gamma=delg*randn(1)+gbar;
%       for j=1:M
%          gam((i-1)*M+j)=gamma;
%          thet((i-1)*M+j)=2*pi*j/M;
%       end
%   end
%  tag = 'uniform in theta, gaussian in gamma';
elseif(iopt==5)
  M=4;  % number of particles in each beamlet
  nb= round(n/M); % number of beamlet via Fawley between 64 to 256 (x16=1024 to 4096)
  if(M*nb~=n)
	error('n must be a multiple of 4')
  end

  % make hammersley pseudo-random sequence of size 1xnb
  dim_num = 1; step = 0; seed = 1; leap = 1; base = 3;
  ham_rand = i_to_hammersley_sequence(dim_num, nb, step, seed, leap, base);
 
  effnoise = sqrt(3*M/(Ns/nb));	% Penman algorithm for Ns/nb >> M
  for i=1:nb
    gamma=delg*randn(1)+gbar;
      for j=1:M
         gam((i-1)*M+j)=gamma;
         thet((i-1)*M+j)=2*pi*j/M + 2*pi*ham_rand(i) + 2*rand(1)*effnoise;
      end
  end

% older loading without Hammersley
%   for i=1:nb
%     gamma=delg*randn(1)+gbar;
%       for j=1:M
%          gam((i-1)*M+j)=gamma;
%          thet((i-1)*M+j)=2*pi*j/M+2*rand(1)*effnoise;
%       end
%   end
%  for j = 1:n
%	gam(j) = gbar+delg*randn(1);
%	thet(j) = 2*pi*j/n + 2*rand(1)*effnoise;
%  end
%  tag = 'Shot Noise in theta, gaussian in gamma';
end

% for hammersley, you need to match the lengths of seed, leap, and base to the dim you
% choose.  For example, try:
% 
% dim_num = 2; n = 100; step = 0; seed = [1,1]; leap = [1,1]; base = [2,3];
% vect = i_to_hammersley_sequence(dim_num,n,step,seed,leap,base);
% vectr = rand(2,n);
% 
% figure(1),plot(vect(1,:),vect(2,:),'.')
% figure(2),plot(vectr(1,:),vectr(2,:),'.')
% Daniel says base=[-n,3] is even better