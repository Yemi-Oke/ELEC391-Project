%%  Variable Declarations

Km=4.44*10^-3; %motor constant (constant) Nm/A [AmberValue: 35.9*10^-3]
nom_current=0.102; %nominal current (data sheet) -> A
nom_speed=1100*pi/30; %nominal speed (data sheet) -> rad/s
J=2.18*10^-7; %rotor inertia (data sheet) -> Nm???
B=7*10^-3;%Km*nom_speed*nom_current; % Mini-Proiject Value-> Nm [7*10^-3] no load torque no load speed
L= 0.026*10^-3; %induction (data sheet) -> henrys [AmberValue: 0.881*10^-3]
R= 0.43; %resitance (data sheet) -> ohms [AmberValue: 16.3]

%SCARA Arm Specifications
l1=98.88*10^-3; %distance from joint1 to joint2 -> m
l2=67.14*10^-3; %distance from joint2 to joint3 -> m

%Control Frequency (already dooubled)
f=25.04; %Hz (see C-Code for derivation)

%DAQ Information
count_rate=1/f; %(from data-sheet)-> this is a placeholder value

%% Deriving the Amplifier Model Transfer Function

s=tf('s');

prompt="1 - MATLAB Model, 2 - SimulationX Mode: ";
rsp=input(prompt);
if(rsp==2)
    G_em=Km/(s*(L*s+R)); 
else
    G_em=Km/(s*(L*s+R)*(J*s+B)); %model of electro-mechanical component
end

H=f/(s+f);%model of the sensor filter; gains are unity so i can ignore them
T=G_em/(1+G_em*H); %system transfer function w/out Amplifier or PID Controller

x=stepinfo(T); %this gives me the step response information -> can also use step(T) for visualization
OS=(x.Peak-x.SettlingMin); % the overshoot=(peak - SS value)
KDC=x.SettlingMin; %the steady state value is the DC gain
pctos=(OS/KDC); %percent overshoot (decimal value)
zeta=sqrt((log(pctos/KDC)^2)/(pi^2+(log(pctos/KDC)^2)));
beta=sqrt(1-zeta^2);

wn=(1/(beta*x.RiseTime))*(pi-atan(beta/zeta)); %Natural Frequency was calculated using the Rise Time

%{
%% Natural Frequency Testing

%step(T);
prompt="1-Setttling Time, 2- Peak Time, 3- Rise Time: "; %comparison use only, will remove later. Just used to see which approximation gives the best Amplifier
choice=input(prompt);
if(choice==2)
    wn=pi/(x.PeakTime*beta); %very noisy
elseif (choice==3)
    wn=(1/(beta*x.RiseTime))*(pi-atan(beta/zeta));
else
    wn=4/((x.SettlingTime)*(zeta));%default to settling time -> very noisy
end
%}

AmplifierTransferFunction=KDC*wn^2/(s^2+2*zeta*wn*s+wn^2);

%[AmberValue: AmplifierTransferFunction=5000/s+5000];
%% Heuristic Tuning & PID Controller Parameters:
p=f; %location of controller pole
G=G_em*AmplifierTransferFunction; %Forward/Open Loop Transfer Function
[Gm_start,Pm_start,Wcg_start,Wcp_start]=margin(G*H); %Wcg represents the XOVER frequency

z=Wcg_start/10; %location of initial controller zero
z_old=0;
Ki_test=1; %integration constant

index=0;
while index<20 %iterate until z value stabilizes
    Kp_test=2/z-2/p;
    Kd_test=1/z^2-Kp_test/p;
    D=Ki_test/s+Kp_test+Kd_test*(s/(s+f)); %PID Contoller Transfer Function
    [Gm_test, Pm_test, Wcg_test, Wcp_test]=margin(D*G*H); %returns Gain Margin, Phase Margin, gain crossover, phase crossover, respectively
    z=Wcg_test/10; 
    index=index+1;
end

%Root Locus
rlocus(D*G*H)
K=1/abs(freqresp(D*G*H,z)); %initial gain K  -> need to modify K value

%{
%Nyquist Plots
nyqlog(D*G*H)
nyqlog(K*D*G*H)
%}

%Closed Loop Transfer Function
CLTF=(K*D*G)/(1+K*D*G*H);
c_step=stepinfo(CLTF);
Kf=c_step.SettlingMin; 

%Final Controller Parameters -> shoud all the constants should be a
%function of Kf
Ki=0;%Kf*Ki_test; %Balance overshoot & steady-state error 
Kp=1;%Kf*Kp_test; %Balance rise time & stability 
Kd=0; %Kf*Kd_test; %Maximize stability 

