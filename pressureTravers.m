clc
clear
close all

fprintf('in the name of God\nHi my name is Behnam and wellcom to my program...\nPlease follow the steps below...\n')

%% inputs ...

L(1)=input('please enter L1 :');
P(1)=input('please enter P1:');
Pdiff=input('set pressure diffrence:');
eLd=input('please estimate length deffrence:');
TL=input('please enter total length:');
RP=input('please enter Reservoir Pressure (Psig):');
OilR=input('please enter Oil Rate (STB/day):');
InsD=input('please enter Instant Well Depth (ft):');
WHP=input('please enter Wellhead Pressure (Psig);');
PD=input('please enter Estimated Pressure Drop (Psig)');
Angle=input('please enter Pipe Inclination Angle with Horizontals:');

limit=input('please enter limit in case you have entered wrong data:');
err='false';
done1='false';
done2='false';
i=1;
%big loop start!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
while strcmp(done1,'false')
    ITER=0;
    AVP=P(i)+PD/2;
    %small loop start***********************************
        while strcmp(done2,'false')
           prop=property(RP,OilR,InsD,WHP,PD,eLd,Angle);
           psessureDropT=prop(1) ;
           PD=prop(2);
           AVT=prop(3);
           cLd=PD/psessureDropT;
           
           if (abs(eLd-cLd)<10^-3)
               break;
           else
              % check err start 
               if ITER>limit
                   err='true';
               else
                   ITER=ITER+1;
                   eLd=cLd;
               end
               % check err end
           end
           
           
        end
    %small loop end************************************
    if strcmp(err,'true')
        disp('somthing went wrong , mybe your data are wrong ,check again...')
        break;
    end
    
    L(i+1)=L(i)+eLd;
    P(i+1)=P(i)+PD;
    
    if L(i+1)<TL
        i=i+1;
    else
        % interpolate 
        p=(P(i+1)+P(i))/2;
        l=(L(i+1)*p-L(i+1)*P(i)-L(i)*p+L(i)*P(i)+L(i)*P(i+1)-L(i)*P(i))/(P(i+1)-P(i));
        fprintf('results:\nl=%s\np=%s\neLd=%s\nPD=%s\n',l,p,eLd,PD);
        break;
    end
    
end
%big loop end!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
