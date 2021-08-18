function [ ans ]=property(RP,OilR,InsD,WHP,PD,eLd,Angle)
% Pressure Drop Calculation
% Summary of this function goes here
% To calculate the pressure drop along the production tubing by using Beggs andBrill correlation
% Detailed explanation goes here
% Input parameters in this function:
% IRP Initial Reservoir Pressure (Psig)
% RP Reservoir Pressure (Psig)
% OilR Oil Rate (STB/day)
% InsD Instant Well Depth (ft)
% WHP Wellhead Pressure (Psig)
% PD Estimated Pressure Drop (Psig)
% Output parameters in this function:
% PDTrue True Pressure Drop (Psig)
% eLd estimate length diffrence
% Angle Pipe Inclination Angle with Horizontal

% Basic Input Parameters

IRT = 180; %Intial Reservoir Temperature (F)
WHT = 120; %Wellhead Temperature (F)
BHT = 180; %Bottomhole Temperature (F)
TD = 8200; %Total Depth (ft)
Dia = 5.5; %Inner Diameter (in)
Rough = 0.0006; %Pipe Roughness (in)
GasG = 0.62; %Gas Gravity
OilAPI = 38.5; %Oil API
OilG = 141.5/(131.5+OilAPI);%Oil Gravit

% Basic Input Parameters End
SD = eLd;% Step Distance (ft)
BHP =WHP; %Bottomhole Pressure (Psig)
AVP = (WHP+BHP)/2;%Average Pressure (Psig)
AVT = ((WHT + (BHT-WHT)/TD*InsD)+(WHT + (BHT-WHT)/TD*(InsD+SD)))/2;%Average Temperature (F)

IIGOR = GasG*(((RP+14.7)/18.2 + 1.4)* PD+10^(0.0125* OilAPI-0.00091*IRT))^(1/0.83);%Instant Intial GOR (scf/STB)

if GasG*(((AVP+14.7)/18.2 + 1.4)*10^(0.0125* OilAPI-0.00091*AVT))^(1/0.83)>IIGOR
    GOR = IIGOR;% Solution GOR (scf/STB)
else
    GOR = GasG*(((AVP+14.7)/18.2 + 1.4)*10^(0.0125* OilAPI-0.00091*AVT))^(1/0.83);% Solution GOR (scf/STB)
end


if GOR<IIGOR
    Bo = 0.972+1.47*10^(-4)*(1.25*AVT+GOR*(GasG/OilG)^0.5)^1.175;% Oil Formation Volume Factor (bbl/STB)
else
    OilC = 10^(-5)*(-1433+5*IIGOR+17.2*AVT-1180*GasG+12.61*OilAPI)/(RP+14.7); %Oil Compressibility
    Bo = (0.972+1.47*10^(-4)*(1.25*AVT+GOR*(GasG/OilG)^0.5)^1.175)*exp(OilC*(RP-AVP));% Oil Formation Volume Factor (bbl/STB)
end


%Z factor and Gas Formation Volume Factor Calculation
Ppc = 756.8-131.07*GasG-3.6*GasG^2;%Pseudocritical Pressure
Tpc = 169.2+349.5*GasG-74*GasG^2; %Pseudocritical Temperature
Pr = (AVP+14.7)/Ppc; %Pseudoreduced Pressure
Tr = (AVT+460)/Tpc; %Pseudoreduced Temperature
A = 1.39*(Tr-0.92)^0.5-0.36*Tr-0.101;
B =(0.62-0.23*Tr)*Pr+(0.066/(Tr-0.86)-0.037)*Pr^2+0.32/10^(9*(Tr-1))*Pr^6;
C = 0.132-0.32*log10(Tr);
D = 10^(0.3106-0.49*Tr+0.1824*(Tr^2));
Z = A+(1-A)/exp(B)+C*Pr^D;% Z Factor
Bg = 0.0282793*(AVT+460)*Z/(AVP+14.7);% Gas Formation Volume Factor (rcf/scf)
GasD = 28.97*GasG*(AVP+14.7)/Z/10.73/(460+AVT);% Gas Density (lb/ft3)


if GasG*(((AVP+14.7)/18.2 + 1.4)*10^(0.0125* OilAPI-0.00091*AVT))^(1/0.83)>IIGOR
    OilC = 10^(-5)*(-1433+5*IIGOR+17.2*AVT-1180*GasG+12.61*OilAPI)/(RP+14.7);
    OilD = (350*OilG+0.0764*GasG*GOR)/5.615/Bo*exp(OilC*(AVP-RP));% Oil Density (lb/ft3)
else
    OilD = (350*OilG+0.0764*GasG*GOR)/5.615/Bo;% Oil Density (lb/ft3)
end


LiqD = OilD;%Liquid Density (lb/ft3)
GasFR = (IIGOR-GOR)* OilR * Bg/86400;% Gas Flowrate (ft3/s)
LiqFR = OilR*Bo*5.615/86400;% Liquid Flowrate (ft3/s)
Area = pi*(Dia/2*0.0833333)^2;% Pipe Sectional Area (ft2)
SfGasFR = GasFR/Area;% Superficial Gas Flowrate (ft/s)
SfLiqFR = LiqFR/Area;% Superficial Liquid Flowrate (ft/s)
SfMixFR = SfGasFR + SfLiqFR;% Superficial Mix Flowrate (ft/s)
GasWFR = GasD * SfGasFR;% Gas Weight Flux Rate (lb/ (ft2.s))
LiqWFR = LiqD * SfLiqFR;% Liquid Weight Flux Rate (lb/ (ft2.s))
MixWFR = GasWFR + LiqWFR;% Mix Weight Flux Rate (lb/ (ft2.s))
NoSlipHUp = LiqFR/(LiqFR + GasFR);% No-slip Holdup
Nfr = SfMixFR^2/32.174/(Dia*0.0833333);% Froude Number


%Gas Viscosity
Mg = 28.967 * GasG;% Mole weight
K1 = (0.00094+2*10^(-6)*Mg)*(AVT+460)^1.5/(209+19*Mg+(AVT+460));
X = 3.5+986/(460+AVT)+0.01*Mg;
Y = 2.4-0.2*X;
GasVisc = K1*exp(X*(GasD/62.4)^Y); %Gas Viscosity (cp)
DeadOilVisc = 10^(AVT^-1.163*exp(6.9824-0.04658*OilAPI))-1;% Dead Oil Viscosity (cp)
A1 = 10.715*(GOR+100)^-0.515;
B1 = 5.44*(GOR+150)^-0.338;


if GasG*(((AVP+14.7)/18.2 + 1.4)*10^(0.0125* OilAPI-0.00091*AVT))^(1/0.83)>IIGOR
    C1 =2.6*AVP^1.187*exp(-11.513+(-8.98*10^-5*AVP));
    OilVisc = A1* DeadOilVisc^B1*(AVP/RP)^C1; % Bubblepoint Oil Viscosity (cp)
else
    OilVisc = A1* DeadOilVisc^B1; % Bubblepoint Oil Viscosity (cp)
end

LiqVisc = OilVisc;
OilIntT = (1.17013-1.694*10^(-3)*AVT)*(38.085-0.259*OilAPI);% Oil Interfacial Tension (dynes/cm)
LiqIntT = OilIntT;% Liquid Interfacial Tension (dynes/cm)
MixVisc = 6.72*10^(-4)*(LiqVisc*NoSlipHUp + GasVisc*(1-NoSlipHUp));% Mix Viscosity (lb/ft/s)
Nre = MixWFR*Dia*0.0833333/MixVisc;% No Slip Reynold Number
Nlv = 1.938*SfLiqFR*(LiqD/LiqIntT)^0.25;% Liquid Velocity Number


%Determine the flow patten
L1 = 316*NoSlipHUp^0.302;
L2 = 0.0009252*NoSlipHUp^-2.4684;
L3 = 0.1*NoSlipHUp^-1.4516;
L4 = 0.5*NoSlipHUp^-6.738;


if (NoSlipHUp<0.01 && Nfr<L1) || (NoSlipHUp>=0.01 && Nfr<L2)
    Flow = 'Segregated';
elseif (NoSlipHUp>=0.01 && NoSlipHUp<0.4 && Nfr>L3 && Nfr<=L1) || (NoSlipHUp>=0.4 && Nfr>L3 && Nfr<=L4)
    Flow = 'Intermittent';
elseif (NoSlipHUp<0.4 && Nfr>=L4) || (NoSlipHUp>=0.4 && Nfr>L4)
    Flow = 'Distributed';
elseif (Nfr>L2 && Nfr<L3)
    Flow = 'Transition';
else
    Flow = ' ';
end

%Liquid Holdup for Horizontal Flow
SegHup0 = 0.98*NoSlipHUp^0.4846/Nfr^0.0868; % Segregated
InterHup0 = 0.845*NoSlipHUp^0.5351/Nfr^0.0173; % Intermittent
DisHup0 = 1.065*NoSlipHUp^0.5824/Nfr^0.0609; % Distributed

% Selected Liquid Holdup
if strcmp(Flow,'Segregated')
    Sel_Hup0 = SegHup0;
elseif strcmp(Flow,'Intermittent')
    Sel_Hup0 = InterHup0;
elseif strcmp(Flow,'Distributed')
    Sel_Hup0 = DisHup0;
else
    Sel_Hup0 = 1;
end

% Final Selected Liquid Holdup
if Sel_Hup0 > NoSlipHUp
    FSel_Hup0 = Sel_Hup0;
else
    FSel_Hup0 = NoSlipHUp;
end

% Correction Factor Coefficient
SegCF = (1-NoSlipHUp)*log(0.011*NoSlipHUp^(-3.768)*Nlv^3.539*Nfr^(-1.614)); % Segregated Flow
InterCF = (1-NoSlipHUp)*log(2.96*NoSlipHUp^(0.305)*Nlv^(-0.4473)*Nfr^(0.0978)); % Intermittent Flow
DisCF = 0; % Distributed Flow


% Selected Correction Factor Coefficient
if strcmp(Flow,'Segregated')
    Sel_CF = SegCF;
elseif strcmp(Flow,'Intermittent')
    Sel_CF = InterCF;
else
    Sel_CF = DisCF;
end

% Final Selected Correction Factor Coefficient
if Sel_CF <= 0
    FSel_CF = 0;
    
else
    FSel_CF = Sel_CF;
end

Rad_Angle = Angle*pi/180; % Pipe Inclination Angle with Horizontal
B_Angle = 1+FSel_CF*(sin(1.8*Rad_Angle)-(1/3)*sin(1.8*Rad_Angle)^3); % Correction Factor

%Liquid Holdup for Transition Flow
AA = (L3-Nfr)/(L3-L2);
BB = 1-AA;

SegB_Angle = 1+SegCF*(sin(1.8*Rad_Angle)-(1/3)*sin(1.8*Rad_Angle)^3); % Correction Factor For Segregated Flow
InterB_Angle = 1+InterCF*(sin(1.8*Rad_Angle)-(1/3)*sin(1.8*Rad_Angle)^3); %Correction Factor For Intermittent Flow

if SegB_Angle<=1
    SegHup1 = SegHup0;
else
    SegHup1 = SegHup0 * SegB_Angle;
end
if InterB_Angle<=1
    InterHup1 = InterHup0;
else
    InterHup1 = InterHup0 * InterB_Angle;
end
TransHup1 = AA* SegHup1 + BB* InterHup1;

% Selected Liquid Holdup for Transition Flow
if TransHup1<=NoSlipHUp
    Sel_TransHup1 = NoSlipHUp;
else
    Sel_TransHup1 = TransHup1;
end

% Final Liquid Holdup For Next Calculation
if strcmp(Flow,'Transition')
    Final_LiqHup = Sel_TransHup1;
elseif B_Angle<=1
    Final_LiqHup = FSel_Hup0;
else
    Final_LiqHup = FSel_Hup0 * B_Angle;
end


TwoPhaseD = Final_LiqHup * LiqD + (1- Final_LiqHup)*GasD; % two-phase density
PressDropG = TwoPhaseD * sin(Rad_Angle)/144; % hydrostatic pressure gradient
% Pressure drop due to friction loss
% Ratio of two-phase to no-slip friction factor calculation

Y = NoSlipHUp/Final_LiqHup^2;
In_Y = log(Y);
S = In_Y/(-0.0523+3.182*In_Y-0.8725*In_Y^2+0.01853*In_Y^4);
if Y == 0
    S_Final = 0;
elseif Y>1 && Y<1.2
    S_Final = log(2.2*Y-1.2);
else
    S_Final = S;
end

FrictionFR = exp(S_Final); % Ratio of two-phase to no-slip friction factor

% No-slip friction factor calculation/ Fanning friction factor calculation
Epsilon_D = Rough/Dia;
FF1 = 64/Nre;
FF2 = 1/(-2*log10(Epsilon_D/3.7+2.51/(Nre*sqrt(FF1))))^2;
FF3 = 1/(-2*log10(Epsilon_D/3.7+2.51/(Nre*sqrt(FF2))))^2;
FF4 = 1/(-2*log10(Epsilon_D/3.7+2.51/(Nre*sqrt(FF3))))^2;
FF5 = 1/(-2*log10(Epsilon_D/3.7+2.51/(Nre*sqrt(FF4))))^2;
FF6 = 1/(-2*log10(Epsilon_D/3.7+2.51/(Nre*sqrt(FF5))))^2;
FF7 = 1/(-2*log10(Epsilon_D/3.7+2.51/(Nre*sqrt(FF6))))^2;
FF8 = 1/(-2*log10(Epsilon_D/3.7+2.51/(Nre*sqrt(FF7))))^2;
FF9 = 1/(-2*log10(Epsilon_D/3.7+2.51/(Nre*sqrt(FF8))))^2;
FF10 = 1/(-2*log10(Epsilon_D/3.7+2.51/(Nre*sqrt(FF9))))^2;
FF11 = 1/(-2*log10(Epsilon_D/3.7+2.51/(Nre*sqrt(FF10))))^2;
FF12 = 1/(-2*log10(Epsilon_D/3.7+2.51/(Nre*sqrt(FF11))))^2;
FF13 = 1/(-2*log10(Epsilon_D/3.7+2.51/(Nre*sqrt(FF12))))^2;
FF14 = 1/(-2*log10(Epsilon_D/3.7+2.51/(Nre*sqrt(FF13))))^2;
FF15 = 1/(-2*log10(Epsilon_D/3.7+2.51/(Nre*sqrt(FF14))))^2;

if Nre<=2100
    FannFriF = FF1 /4; % Fanning friction factor
else
    FannFriF = FF15 /4;
end

TwoPhaseFriF = FrictionFR * FannFriF; % Two-phase friction factor
PressDropF = 2*TwoPhaseFriF*SfMixFR*MixWFR/144/32.174/(Dia*0.0833333); % pressure gradient due to the friction pressure loss
PressDropT = PressDropG + PressDropF; % Total pressure gradient
PDTrue = SD*PressDropT; % Ture pressure drop in the length increment
ans=[PressDropT,PDTrue,AVT];
end