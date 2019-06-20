function [T0,Angle]=CrystalT0(Reflection,Energy,thickness,Bragg,Eta)
%Reflection [h,k,l], as line vector
%thickness in micron
%Energy in eV
%T0 is output in micron
%Angle in degrees
%Bragg accepts values 'b'=bragg, 'l'=laue
%Eta asymmetry angle

LUT=[1,1,1,3.01034,1.09,8.17,192.0;
    2,2,0,4.9151,1.98,3.04,106.0;
    3,1,1,5.76401,3.74,2.20,56.0;
    4,0,0,6.95161,3.63,1.51,60.6
    7,3,3,14.2251,19.5,0.36,12.6
];

c=299792458;

LUT(:,4)=LUT(:,4)*1000;

[LN,LM]=size(LUT);

for II=1:LN
    if(LUT(II,1:3)==Reflection)
        break
    end
end

Angle=asin(LUT(II,4)/Energy);
Lambdas=LUT(II,5);

if Bragg=='b'
    if(isnan(Eta))
        Eta=0; %Asimmetry Angle Eta=0 Simmetric Bragg Case
    end
    psi0=Angle-pi/2+Eta;%+Eta-pi/2;
    gamma0=cos(psi0);
    sin(Angle);
else
    if(isnan(Eta))
        Eta=pi/2; %Asimmetry Angle Eta=pi/2 Simmetric Laue Case
    end
    psi0=Angle-pi/2+Eta;%+Eta-pi/2;
    gamma0=cos(psi0);
end

effective_length=thickness*10^-6/gamma0;

T0= 2*(Lambdas*10^-6)^2/(c*effective_length);

%2*(Lambdas*10^-6)^2/(c*(thickness*10^-6))

Angle=Angle*360/2/pi;