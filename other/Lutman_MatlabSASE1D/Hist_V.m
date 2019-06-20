function Isto=Hist_V(V,EV,BINS,GA)

MV=mean(V);
if(nargin<4)
    GA=0;
end
NV=V/MV;
Min=min(NV);
Max=max(NV);

if(EV)
    Nodes=linspace(Min,Max,BINS-1);
    Nodes=[Nodes(1)-abs(Nodes(2)-Nodes(1)),Nodes,Nodes(BINS-1)+abs(Nodes(2)-Nodes(1))];
else
    Nodes=logspace(0,log10(Max-Min+1),BINS-1)+Min-1;
    Nodes=[Nodes(1)-abs(Nodes(2)-Nodes(1)),Nodes,Nodes(BINS-1)+abs(Nodes(2)-Nodes(1))];
end

Lost=0;
Events=zeros(1,length(Nodes));
if(GA)
    Resorting=V*0;
end

for JJ=1:length(NV)
   Indice=find(NV(JJ)>Nodes,1,'Last');
   if(Indice>length(Events))
       Lost=Lost+1;
   else
       Events(Indice)=Events(Indice) + 1;
   end
   if(GA)
       Resorting(JJ)=Indice;
   end
end

Lost

MULTI=abs(diff(Nodes));
EventsC=Events(1:(length(Events)-1))./MULTI;
Asse=(Nodes(1:(length(Nodes)-1))+Nodes(2:(length(Nodes))))/2;
Integrale=Trapez(Asse,EventsC)
EventsC=EventsC./Integrale;
media=Trapez(Asse,Asse.*EventsC);

Isto.media=media;
Isto.Asse=Asse;
Isto.EventsC=EventsC;
Isto.Events=Events;
if(GA)
    Isto.Resorting=Resorting;
end