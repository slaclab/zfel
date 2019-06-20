%Look at undulator configuration
SCANAFTERMODULE=str2num(SCANAFTERMODULE);
if(strcmp(ManualUndConf,'No Manual Configuration Stored')) %standard undulators
    if(strcmp(CrystalInString,'No Crystal Inserted'))
        if(BRANCHINGCYCLE==0)
            Element{1}.type='U';
            Element{1}.length=T3{2}*T3{3};
            Element{1}.Kvalue=T3{4}*ones(1,Element{1}.length*T1{1});
        else % a scan is a break point
            Element{1}.type='U';
            Element{1}.length=T3{2}*SCANAFTERMODULE;
            Element{1}.Kvalue=T3{4}*ones(1,Element{1}.length*T1{1});
            Element{2}.type='S';
            Element{2}.length=0;
            Element{2}.Kvalue=T3{4};
            Element{3}.type='U';
            Element{3}.length=T3{2}*(T3{3}-SCANAFTERMODULE);
            Element{3}.Kvalue=T3{4}*ones(1,Element{3}.length*T1{1});
        end
    else % Crystal is in, that is a break point!
        if(BRANCHINGCYCLE==0)
            Element{1}.type='U';
            Element{1}.length=T3{2}*T4{1};
            Element{1}.Kvalue=T3{4}*ones(1,Element{1}.length*T1{1});
            Element{2}.type='C';
            Element{2}.length=0;
            Element{2}.Kvalue=T3{4};
            Element{3}.type='U';
            Element{3}.length=T3{2}*(T3{3}-T4{1});
            Element{3}.Kvalue=T3{4}*ones(1,Element{3}.length*T1{1});
        else %lots of breakpoints!!
            if(SCANAFTERMODULE<T4{1})
                Element{1}.type='U';
                Element{1}.length=T3{2}*SCANAFTERMODULE;
                Element{1}.Kvalue=T3{4}*ones(1,Element{1}.length*T1{1});
                Element{2}.type='S';
                Element{2}.length=0;
                Element{2}.Kvalue=T3{4};
                Element{3}.type='U';
                Element{3}.length=T3{2}*(T4{1}-SCANAFTERMODULE);
                Element{3}.Kvalue=T3{4}*ones(1,Element{3}.length*T1{1});
                Element{4}.type='C';
                Element{4}.length=0;
                Element{4}.Kvalue=T3{4};
                Element{5}.type='U';
                Element{5}.length=T3{2}*(T3{3}-T4{1});
                Element{5}.Kvalue=T3{4}*ones(1,Element{5}.length*T1{1});
            elseif(SCANAFTERMODULE==T4{1})
                Element{1}.type='U';
                Element{1}.length=T3{2}*SCANAFTERMODULE;
                Element{1}.Kvalue=T3{4}*ones(1,Element{1}.length*T1{1});
                Element{2}.type='S';
                Element{2}.length=0;
                Element{2}.Kvalue=T3{4};
                Element{3}.type='C';
                Element{3}.length=0;
                Element{3}.Kvalue=T3{4};
                Element{4}.type='U';
                Element{4}.length=T3{2}*(T3{3}-SCANAFTERMODULE);
                Element{4}.Kvalue=T3{4}*ones(1,Element{4}.length*T1{1});
            else(SCANAFTERMODULE>T4{1})
                Element{1}.type='U';
                Element{1}.length=T3{2}*T4{1};
                Element{1}.Kvalue=T3{4}*ones(1,Element{1}.length*T1{1});
                Element{2}.type='C';
                Element{2}.length=0;
                Element{2}.Kvalue=T3{4};
                Element{3}.type='U';
                Element{3}.length=T3{2}*(SCANAFTERMODULE-T4{1});
                Element{3}.Kvalue=T3{4}*ones(1,Element{3}.length*T1{1});
                Element{4}.type='S';
                Element{4}.length=0;
                Element{4}.Kvalue=T3{4};
                Element{5}.type='U';
                Element{5}.length=T3{2}*(T3{3}-SCANAFTERMODULE);
                Element{5}.Kvalue=T3{4}*ones(1,Element{5}.length*T1{1});
            end
        end
    end
else % undulators with breaks % Many events can break the line and has to restart things
    [ElementNumber,Six]=size(UndulatorManualConfiguration);
%     Undulengths(1)=0;
%     Kvalue{ThisBreakCounter}=[];
    ElementCounter=0;
    UNDIN=0;
    for KK=1:ElementNumber %Reads the structure step by step and determines the number of breaks...
        if(UndulatorManualConfigurationRows{KK}(1)=='U') %it is an undulator
            UNDIN=UNDIN+1;
            if(ElementCounter)
                if(Element{ElementCounter}.type=='U');
                    disp(['mle ',UndulatorManualConfiguration{KK,1}]);
                Element{ElementCounter}.length=Element{ElementCounter}.length+UndulatorManualConfiguration{KK,4};
                if(~isnan(UndulatorManualConfiguration{KK,1}))
                   ThisTaper=ones(1,UndulatorManualConfiguration{KK,4}*T1{1})*UndulatorManualConfiguration{KK,1};
                else %Cont Taper
                    %disp('Continuos taper')
                   TaperData=str2num(UndulatorManualConfiguration{KK,6});
                   ThisTaper=ones(1,UndulatorManualConfiguration{KK,4}*T1{1}).*TaperData(1) + ones(1,UndulatorManualConfiguration{KK,4}*T1{1})/T1{1}*TaperData(2) +  (ones(1,UndulatorManualConfiguration{KK,4}*T1{1})/T1{1}>=TaperData(3)).*((ones(1,UndulatorManualConfiguration{KK,4}*T1{1})/T1{1}-TaperData(3)).^2)*TaperData(4);
                end
                Element{ElementCounter}.Kvalue=[Element{ElementCounter}.Kvalue,ThisTaper]
                else
                    ElementCounter=ElementCounter+1;
                    Element{ElementCounter}.type='U';
                    Element{ElementCounter}.length=UndulatorManualConfiguration{KK,4};
                    if(~isnan(UndulatorManualConfiguration{KK,1}))
                        Element{ElementCounter}.Kvalue=ones(1,UndulatorManualConfiguration{KK,4}*T1{1})*UndulatorManualConfiguration{KK,1};
                    else %Cont Taper
                        TaperData=str2num(UndulatorManualConfiguration{KK,6});
                        Element{ElementCounter}.Kvalue=ones(1,UndulatorManualConfiguration{KK,4}*T1{1}).*TaperData(1) + ones(1,UndulatorManualConfiguration{KK,4}*T1{1})/T1{1}*TaperData(2) +  (ones(1,UndulatorManualConfiguration{KK,4}*T1{1})/T1{1}>=TaperData(3)).*((ones(1,UndulatorManualConfiguration{KK,4}*T1{1})/T1{1}-TaperData(3)).^2)*TaperData(4);
                    end
                end       
            else
                ElementCounter=ElementCounter+1;
                Element{ElementCounter}.type='U';
                UndulatorManualConfiguration{KK,1}
                Element{ElementCounter}.length=UndulatorManualConfiguration{KK,4};
                if(~isnan(UndulatorManualConfiguration{KK,1}))
                    Element{ElementCounter}.Kvalue=ones(1,UndulatorManualConfiguration{KK,4}*T1{1})*UndulatorManualConfiguration{KK,1};
                else %Cont Taper
%                     save Temp4s
%                     disp('Continuos taper')
                    TaperData=str2num(UndulatorManualConfiguration{KK,6});
                    Element{ElementCounter}.Kvalue=ones(1,UndulatorManualConfiguration{KK,4}*T1{1}).*TaperData(1) + ...
ones(1,UndulatorManualConfiguration{KK,4}*T1{1})/T1{1}*TaperData(2) + ...
((1:(UndulatorManualConfiguration{KK,4}*T1{1}))/T1{1}>=TaperData(3)).*(((1:(UndulatorManualConfiguration{KK,4}*T1{1}))/T1{1}-TaperData(3)).^2)*TaperData(4);
                end
            end
            
            if(UndulatorManualConfiguration{KK,3})
                ElementCounter=ElementCounter+1;
                Element{ElementCounter}.type='K';
                Element{ElementCounter}.length=0;
            end
            UndulatorManualConfiguration{KK,2}
            if((UndulatorManualConfiguration{KK,2}))
                ElementCounter=ElementCounter+1;
                Element{ElementCounter}.type='D';
                Element{ElementCounter}.length=UndulatorManualConfiguration{KK,2};
            end
        else
            ElementCounter=ElementCounter+1;
            Element{ElementCounter}.type='C';
            Element{ElementCounter}.length=0;
        end
        if(UNDIN==SCANAFTERMODULE)
            ElementCounter=ElementCounter+1;
            Element{ElementCounter}.type='S';
            Element{ElementCounter}.length=0;
        end
    end    
end