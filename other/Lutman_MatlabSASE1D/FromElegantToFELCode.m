load toALnohead.dat
BunchBins=1000;
BunchSample=Hist_V(toALnohead(:,1),1,BunchBins,1);
TimeWidth=max(toALnohead(:,1))-min(toALnohead(:,1));
TimeAxis=linspace(0,TimeWidth,BunchBins);
for II=1:BunchBins
   InthisBin=find((BunchSample.Resorting==II));
   if(~isempty(InthisBin))
       MeanEnergy(II)=mean(toALnohead(InthisBin,2));
       StdEnergy(II)=std(toALnohead(InthisBin,2));
       ParticleIn(II)=length(InthisBin);
   else
       MeanEnergy(II)=mean(toALnohead(:,2));
       StdEnergy(II)=0;
       ParticleIn(II)=0;
   end
end