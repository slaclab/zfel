% Field Setup
fieldr=zeros(s_steps,1);
fieldi=zeros(s_steps,1);
% ar = zeros(s_steps,z_steps+1);
% ai = zeros(s_steps,z_steps+1);

if(1) %use this if to ovveride this in case of external field
    input_field = sqrt(T5{1})*exp(-((1:T1{2})-T1{2}/2).^2/(2*(5*T1{2}/2)^2));  % scaled initial input field vector [sqrt(W)] 
end

if(SeedShape==1)
    fieldr(:,1) = real(input_field)'/sqrt(rhoPbeam);
    fieldi(:,1) = imag(input_field)'/sqrt(rhoPbeam);
elseif(SeedShape==2) % load seed setup and overrides other paramenters (Amplitude + Phase)
    fieldr(:,1) = max(input_field)
else %use this for external field from file
    fieldr(:,1) = real(input_field)'/sqrt(rhoPbeam);
    fieldi(:,1) = imag(input_field)'/sqrt(rhoPbeam);
end

Final_field=fieldr+1i*fieldi;