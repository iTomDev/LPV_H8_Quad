
% Simple quadrotor model linearised at one operating point, non-PV
% used for synthesising hinf LSDP controllers
% Based on each axis being independent, small angles, 2nd order DEs

% It's not that good but its in the right ballpark for a quadrotor

% Thomas Pile, , SHU
% 21/6/2018

% Jxx = Jyy = 0.25*m*r^2 = 0.25*0.5*0.15^2 = 0.0028
% Jzz = 0.5*m*r^2 = 0.5*0.5*0.02^2 = 1.0000e-04 = 0.0001

% Constants
m = 0.5;
kr = -0.02;
Ixx = 0.0028;
Iyy = 0.0028;
Izz = 0.00110;

% Model
% ddrpy: model of the angular acceleration
% Could be modelled as 3 second order DE
a44 = kr*(1/Ixx);
a55 = kr*(1/Iyy);
a66 = kr*(1/Izz);
% off diagonal
a45 = -(Izz-Iyy)/Ixx;
a46 = -(Izz-Iyy)/Ixx;
a54 = -(Ixx-Izz)/Iyy; 
a56 = -(Ixx-Izz)/Iyy;
a64 = -(Iyy-Ixx)/Izz;
a65 = -(Iyy-Ixx)/Izz;
% d2
ddrpy = [a44 a45 a46;
         a54 a55 a56;
         a64 a65 a66];
     
% drpy: model of the angular rate 
% This is often modelled as an integrator, but actually there is a drag
% component to consider 
% first order part of the 3 DE
a14 = 0.99;
a15 = 0.99;
a16 = 0.99;
drpy = [a14 0   0;
        0   a15 0;
        0   0   a16];
% combined
A = [zeros(3,3) drpy;
     zeros(3,3) ddrpy];

% control inputs: roll distrib, pitch distrib, yaw distrib, total thrust
% scale for torque. 1/0.01Nm which is enough to rotate with unit input
% b1 = (1/Ixx)*0.0065; % so an input of 5 would be about right for 0.05nm 
% b2 = (1/Iyy)*0.012;
% b3 = (1/Izz)*0.003;
% scale so that maximum input (255) is 1
b1 = (1/Ixx) *(0.3060); 
b2 = (1/Iyy) *(0.3060);
b3 = (1/Izz) *(0.3060);
B = [0  0  0  0
     0  0  0  0
     0  0  0  0
     b1 0  0  0
     0  b2 0  0
     0  0  b3 0];
% B = [0  0  0  
%      0  0  0  
%      0  0  0  
%      b1 0  0  
%      0  b2 0 
%      0  0  b3];
%  
% select rpy, xyz as outputs
%C = [eye(3)         zeros(3,3)]; % radians. This is good for a max output of 1 ~57 degrees
C = [57.2958*eye(3) zeros(3,3)];  % degrees. Good for short responses 

D = zeros(3,4);

G = ss(A,B,C,D);

% provisional analysis of model
pzplot(G)
step(G,1)
%title('Pitch Coefficient Model Test')
%ylabel('angle [degrees]')
%xlabel('time [seconds]')

%%

% control synth
s = tf('s');
W1 = tf([1.3],[0.5 5.9]);
W2 = 0.4;
sigma(G,W1*G*W2)
legend('G','W1*G*W2')
title('Shaped Plant')

[K,CL,GAM,INFO] = ncfsyn(G,W1,W2)
sigma(K)

K = -K;

step(feedback(G,K),3)
step(CL((4:6),(4:6)),5)
% loop functions
loops = loopsens(G,K);
sigma(loops.So,loops.To,loops.CSo);
legend('S','T','KS');
title('Loop Shapes: S, T and KS')
% controller poles
pzplot(-K)
title('Controller Poles')

% TF from ref to error (I+KG)^-1 = Si
sigma(loops.Si)

% controller order reduction
% reduce controller order
[K2, redinfo2] = reduce(K,9);
[Ka Kb Kc Kd] = ssdata(K)

% export controller
%[Ka Kb Kc Kd] = ssdata(K);
save('Ka.mat','Ka');
save('Kb.mat','Kb');
save('Kc.mat','Kc');
save('Kd.mat','Kd');
size(Ka)
K2 = inv(W1)*K;

% reduced order controller
step(feedback(G,K),'b--',feedback(G,K2),'r--',3)
legend('Full Order','7th Order Controller')
title('Step Response: Reduced Order')

% loop functions
loops_reduced = loopsens(G,K2);
sigma(loops_reduced.So,loops_reduced.To,loops_reduced.CSo);
legend('S','T','KS');
title('Loop Shapes: S, T and KS - Reduced Order')

%%
% Uncertainty modelling
Ixx = ureal('Ixx',0.0028,'percentage',[-10 10]);
Iyy = ureal('Iyy',0.0028,'percentage',[-10 10]);
Izz = ureal('Izz',0.00110,'percentage',[-10 10]);

a44 = kr*(1/Ixx);
a55 = kr*(1/Iyy);
a66 = kr*(1/Izz);
a45 = -(Izz-Iyy)/Ixx;
a46 = -(Izz-Iyy)/Ixx;
a54 = -(Ixx-Izz)/Iyy; 
a56 = -(Ixx-Izz)/Iyy;
a64 = -(Iyy-Ixx)/Izz;
a65 = -(Iyy-Ixx)/Izz;
ddrpy = [a44 0   0;
         0   a55 0;
         0   0   a66];
ddrpy = [a44 a45 a46;
         a54 a55 a56;
         a64 a65 a66];
a14 = 0.99;
a15 = 0.99;
a16 = 0.99;
drpy = [a14 0   0;
        0   a15 0;
        0   0   a16];
A = [zeros(3,3) drpy;
     zeros(3,3) ddrpy];
b1 = (1/Ixx) *(0.3060); 
b2 = (1/Iyy) *(0.3060);
b3 = (1/Izz) *(0.3060);
B = [0  0  0  0
     0  0  0  0
     0  0  0  0
     b1 0  0  0
     0  b2 0  0
     0  0  b3 0];
C = [57.2958*eye(3) zeros(3,3)];  % degrees. Good for short responses 
D = zeros(3,4);
% form uss model
uG = uss(A,B,C,D)
% verify nominal model
CLP = feedback(uG.NominalValue,K);
step(CLP,5)
% uncertain model
CLP = feedback(uG,K);
% uncertain time response
step(usample(CLP,10))
title('Uncertain Step Response')
% uncertain singular values
sigma(usample(CLP,10))
title('Uncertain Singular Value Plot of CL')

% mu analysis
LowerBound: 7.1559
UpperBound: 10.0000
CriticalFrequency: Inf
[stabmarg,wcu] = robstab(CLP); 

