%*************************************************************************
%  
%  (Mulitple) paramter-dependent robust controller design for miniaturized optical image stabilizers 
%             ----------------------------------------
%
% By Pan Zhao, Control Engineering Lab, UBC, BC, Canada
% Under supervision of Dr. Ryozo Nagamune.
% Revision: Feb 25, 2018
% Creation: Jan 16, 2016.

% Function: Main function for designing (single and multiple) parameter-depedent robust controllers
% Reference: Zhao, R. Nagamune and M. Chiao, \Multiple parameter-dependent robust control of miniaturized optical
% image stabilizers," accepted by Control Engineering Practice, 2018.

clear 
% close all;
%% Parameters for controller design
w20 = 101.7*2*pi;   
XY_PD = 3;                          % 0 for PiD X a nd Y, 1 for PD X and constant Y; 2 for PD Y and constant X; 3 for PD X and Y

% wn2 for five plants are: 97.11, 101.6, 101.7,104.9,106.2 
Theta1 = [95 103; 104 108]*2*pi;   % use [95 108]*2*pi for design a single PD robust controler

Theta2Part = 0;% 1 for partition of theta_2
% for uncertainty in w2n ,and removing filter, epsi should be >1e10 
Epsi= logspace(2,3,5);%line search parameter in Sato's method, logarithmically spaced
Epsi = Epsi(2);
% Epsi = 110;          

%% Plant parameters
DCgain = 10^(12.29/20); %12.41;
zeta1 = 5e-3; 
% w10 = 64.379*2*pi;%63.84

zeta2 = 5e-4;
Td = 21.86/(360*52.87); % pure time delay

%% Uncertainty in measurement
deltam = 1*2*pi %0.1 *2*pi
Delta = [-deltam, deltam];

LMI0 = [];

%% plant
syms theta1 theta2 theta3 % theta1: w20, theta2: K1
% theta1 = w20;
% theta2 = 0;

K1 = theta3;
w10 = theta2; 
K2 = DCgain;
% K1 = 0.1;K2 = DCgain; theta1 = w20;
%%  the 1st mode

% A1 = [0 1; -w10^2/w20^2*theta1^2 -2*zeta1*w10/w20*theta1];
% B1 = [0 1]';
% C1 = [K1*w10^2/w20^2*theta1^2 0];
% D1 = 0;

% to make C2 constant
% A1 = [0 K1*w10^2/w20^2*theta1^2; -1/K1 -2*zeta1*w10/w20*theta1];
% B1 = [0 1]';
% C1 = [1 0];
% D1 = 0;

% theta1 = w20;
%  Another way to improve numerical stability
% A1 = [0 w10; -w10/w20^2*theta1^2 -2*zeta1*w10/w20*theta1]; %use 1 instead of /w20*theta1  to make B term constant
% B1 = [0 1]';
% C1 = [K1*w10/w20^2*theta1^2 0];
% D1 = 0;

% % to remove theta1^2 term
% A1 = [0 w10/w20*theta1; -w10/w20*theta1 -2*zeta1*w10/w20*theta1]; %use 1 instead of /w20*theta1 to make B term constant
% B1 = [0 1]';
% C1 = [K1*w10/w20*theta1 0];
% D1 = 0;

% To treat w10 as a seperate GS para
A1 = [0 w10; -w10 -2*zeta1*w10]; %use 1 instead of /w20*theta1 to make B term constant
B1 = [0 1]';
C1 = [K1*w10 0];
D1 = 0;

% to separate K1 and theta1
% A1 = [0 w10/w20^2*theta1^2; -w10 -2*zeta1*w10/w20*theta1]; %use 1 instead of /w20*theta1  to make B term constant
% B1 = [0 1]';
% C1 = [K1*w10 0];
% D1 = 0;

% A1 = [0 w10/w20*theta1; -w10/w20*theta1 -2*zeta1*w10/w20*theta1]; %use 1 instead of /w20*theta1  to make B term constant
% B1 = [0 K1]';
% C1 = [1*w10/w20*theta1 0];
% D1 = 0;

% use theta1/w20 as the parameter 
% A1 = [0 w10; -w10*theta1^2 -2*zeta1*w10*theta1]; %use 1 instead of /w20*theta1  to make B term constant
% B1 = [0 1]';
% C1 = [K1*w10*theta1^2 0];
% D1 = 0;
%% 2nd mode
% A2 = [0 1;-theta1^2 -2*zeta2*theta1];
% B2 = [0 1]';
% C2 = [K2*w20^2 0];
% D2 = 0;

% another way to improve numerical stability
A2 = [0 w20;-theta1^2/w20 -2*zeta2*theta1];% use w20 instead of theta1 to make B term constant
B2 = [0 1]'*40;
C2 = [K2*w20 0]/40;
D2 = 0;

% use theta1/w20 as the parameter 
% A2 = [0 w20;-w20*theta1^2 -2*zeta2*w20*theta1];% use w20 instead of theta1 to make B term constant
% B2 = [0 1]'*100;
% C2 = [K2*w20 0]/100;
% D2 = 0;

% pade approximation of pure delay
% A3 = -2/Td;
% B3 = 1;
% C3 = 4/Td;
% D3 = -1;
tf3 = ss(tf([-0.5*Td  1],[0.5*Td  1]));
A3 = tf3.a; B3 = tf3.b; C3 = tf3.c; D3 = tf3.d;

% sum of (A1,B1,C1,D1) and (A2,B2,C2,D2)
A = [A1 zeros(size(A1,1),size(A2,2)); zeros(size(A2,1),size(A1,2)) A2];
B = [B1;B2];
C = [C1 C2];
D = D1+D2;

% only consider the second mode
% A = A2; B = B2; C = C2; D = D2;

% only consider the first mode
% A = A1; B = B1; C = C1; D = D1;

% product of (A,B,C,D) and (A3,B3,C3,D3)
% A = [A B*C3; zeros(size(A3,1),size(A,2)) A3];
% B = [B*D3; B3];
% C = [C D*C3];
% D = D*D3;

%% check whether it is equal to the transfer function mdoel
% TF_SS = tf(ss(A,B,C,D))
% mode1 = tf(K1*(w10/w20*theta1)^2,[1 2*zeta1*w10/w20*theta1 (w10/w20*theta1)^2]);
% mode2 = tf(K2*w20^2,[1 2*zeta2*theta1 theta1^2]);
% [numD,denD] = pade(Td,1);
% mode3 = tf(numD,denD);
% TF = (mode1+mode2)*mode3
% hold on;
% bode(TF_SS,TF);return;

%% 

% %% another way for designing Hinf controller
% K1 = 0.1;K2 = DCgain;
% theta2 = w20;
% Model1 = K1*tf(w10^2,[1 2*zeta1*w10 w10^2]);
% Model2 = K2*tf(w20^2,[1 2*zeta2*w20 w20^2]);
% DelayT = 21.86/(360*52.87);
% [numD,denD] = pade(DelayT,1);
% DelayMod = tf(numD,denD);
% Model = (Model1+Model2)*DelayMod; 
% sys = ss(Model); 
% A = sys.a;
% B = sys.b;
% C = sys.c;
% D = sys.d;
n = size(A,1);

%% Generalized plant 
Weinv = makeweight(0.028,750,3);%(0.001,20,2); 
We = 1/5*ss(inv(Weinv));
% We = We*ss(tf(1,[1e-5 1]))*ss(tf(1,[1e-5 1]))*ss(tf(1,[1e-5 1]));
ne = order(We);
%% penalize the input around the 2nd mode
% Wu = 1/100*1/5*ss(11.5*tf([1/(w20*0.5) 1],[1/(w20*0.7) 1])*tf([1/(w20*1.4) 1],[1/(w20*1) 1])); 
% nu = order(Wu);
% Gasym.A = [A zeros(n,ne+nu); -We.b*C We.a zeros(ne,nu); zeros(2,n+ne) Wu.a];
% Gasym.B1 = [zeros(n,1); We.b ;zeros(2,1)]; Gasym.B2 = [B; -We.b*D; Wu.b];
% Gasym.C1 = [-We.d*C We.c zeros(1,nu); zeros(1,n+ne) Wu.c];
% Gasym.D11 = [We.d; zeros(1,1)]; Gasym.D12 =[-We.d*D; Wu.d];
% Gasym.C2 = [-C zeros(1,ne+nu)]; Gasym.D21 = 1; Gasym. D22 = -D;

Wu = 0.023; 
Gasym.A = [A zeros(n,ne); -We.b*C We.a];
Gasym.B1 = [zeros(n,1); We.b]; Gasym.B2 = [B; -We.b*D];
Gasym.C1 = [-We.d*C We.c; zeros(1,n+ne)];
Gasym.D11 = [We.d; zeros(1,1)]; Gasym.D12 =[-We.d*D; Wu];
Gasym.C2 = [-C zeros(1,ne)]; Gasym.D21 = 1; Gasym. D22 = -D;
Wu = ss(Wu);
%% Pre filter y to remove parameter dependence of C2
n =  size(Gasym.A,1);
nz = size(Gasym.C1,1);
ny = size(Gasym.C2,1); 
Fy = ss(tf(1,[1/5e3 1]));
Gasym.A = [Gasym.A zeros(n,order(Fy)); Fy.b*Gasym.C2 Fy.a];
Gasym.B1 = [Gasym.B1; Fy.b*Gasym.D21];
Gasym.B2 = [Gasym.B2; Fy.b*Gasym.D22];
Gasym.C1 = [Gasym.C1 zeros(nz,order(Fy))]; 
Gasym.D11 = Gasym.D11; Gasym.D12 = Gasym.D12;
Gasym.C2 = [zeros(ny,n) Fy.c];
Gasym.D21 = 0; Gasym.D22 = 0;

% %% Post filter u to remove parameter dependence of B3
% n =  size(Gasym.A,1); 
% Fu = ss(tf(1,[1/5e3 1])); %Fu.c = Fu.c/10; Fu.b = Fu.b*10;
% Gasym.A = [Gasym.A Gasym.B2*Fu.c; zeros(1,n) Fu.a];
% Gasym.B1 = [Gasym.B1;0]; Gasym.B2 = [zeros(n,1);Fu.b];
% Gasym.C1 = [Gasym.C1 Gasym.D12*Fu.c];
% Gasym.D11 = Gasym.D11; Gasym.D12 = 0;
% Gasym.C2 = [Gasym.C2 0];
% Gasym.D21 = Gasym.D21; Gasym.D22 = 0;

% H inf controller design
% Ga = AugPltEv(Gasym, [w20 0.1]');
% % Ga = Gasym;
% Ga = ss(Ga.A, [Ga.B1 Ga.B2], [Ga.C1; Ga.C2], [Ga.D11 Ga.D12;Ga.D21 Ga.D22]);
% Ga = balred(Ga,order(Ga));
% [Khinf,clp,Gam] = hinfsyn(Ga,1,1,'method','lmi'); Gam

B2 = Gasym.B2; C2 = Gasym.C2;
GSParaNum = 3;
IsGridding = [1 1 0]; % Donnot need to grid the value set of 2nd GS para. since LMI is affine w.r.t it. 

% Fcn_theta = @(x) [1 x(1) x(1)^2]; %function for PD matrices
% d_Fcn_theta = @(x) [0 1 2*x(1)];  %function for derivative of PD matrices
% FthetaNum = [1 2]; %square form    

Fcn_theta = @(x) [1 x(1) x(2)]; %function for PD matrices
d_Fcn_theta = @(x) [0 1 1];  %function for derivative of PD matrices
FthetaNum = [1 1 1]; %affine form  

% Fcn_theta = @(x) [1 x(1) x(2)]; %function for PD matrices
% d_Fcn_theta = @(x) [0 1 1];  %function for derivative of PD matrices
% FthetaNum = [1 1 1]; %affine form 


% Theta2 = [60 70]*2*pi; % for w10 
Theta2 = [61 65; 66 71]*2*pi
Theta3 = [-0.18 0.18];
% Theta3 = 0;
d_thetah = 0;  
theta_min = [Theta1(1,1) Theta2(1,1) Theta3(1,1)]';theta_max = [Theta1(end,2) Theta2(end,2) Theta2(end,2)]';    
NumGain = 1; %for numerical issue
n =  size(Gasym.A,1);
nw = size(Gasym.B1,2);
nu = size(Gasym.B2,2);
nz = size(Gasym.C1,1);
ny = size(Gasym.C2,1); 
%% scalar parameters initialization
if deltam == 0 
    Epsi = 0
else
%      Epsi = Epsi(9);
end 
d_Thetah = {[-d_thetah d_thetah],[-d_thetah d_thetah],[-d_thetah d_thetah]};  
%% regnum and ssnum determination,
regnum1 = size(Theta1,1); %num. of subsets for theta1, only partition theta_1
regnum2 = size(Theta2,1); % num. of subsets for theta2;
regnum = regnum1 * regnum2; % Total num. of subsets. S shape from the bottom to order them, an example:
% 4 5 6
% 1 2 3
SSNum = (regnum1-1)*2*regnum2+(regnum2-1)*2*regnum1; % switching surface num
%% Determination of regid1 & regid2  
REGID = zeros(regnum,3);
for regid = 1:regnum
   if mod(regid,regnum1) == 0  
        regid1 = regnum1;
    else
        regid1 = mod(regid,regnum1);
   end
   if regid > regnum1
        regid2 = 1+ floor(regid/(regnum1+0.1));
    else
        regid2 = 1;
   end  
   REGID(regid,:) = [regid regid1 regid2];   
end
clear regid1 regid2;
%% Paras 
I = eye(n);
%% gridding acquistion
% deltam = 0;
delta_l= -deltam;delta_u = deltam;
for regid1 = 1:regnum1
    [ThetaT{regid1},DeltaT{regid1}] = AdmRegGrid(Theta1(regid1,:),delta_l,delta_u,theta_min(1),theta_max(1),1,0); 
end
delta_l = 0; delta_u = 0;
if GSParaNum > 1
    for regid2 = 1:regnum2
        [Theta2T{regid2},Delta2T{regid2}] = AdmRegGrid(Theta2(regid2,:),delta_l,delta_u,theta_min(2),theta_max(2),1,0);
        if sum(abs(Delta2T{regid2})) == 0
            Theta2T{regid2} = unique(Theta2T{regid2});
            Delta2T{regid2} = Theta2T{regid2}*0;
        end
    end
end

   
% % Theta2T{1} (Theta2T{1}==0) = 1e-3;
Theta3T{1} = Theta3;
Delta3T{1} = 0-Theta3T{1}; % note that Delta2T is not determined by deltam   
% Theta2T{1} (Theta2T{1}==0) = 1e-3;
% Theta2T{1} = 0;
% Delta2T{1} = 0; % note that Delta2T is not determined by deltam 

TD.ThetaT = ThetaT;
TD.DeltaT = DeltaT;
TD.Theta2T = Theta2T; 
TD.Delta2T = Delta2T;
TD.Theta3T = Theta3T; 
TD.Delta3T = Delta3T;
TD.Theta1 = Theta1;
TD.Theta2 = Theta2;
TD.Theta3 = Theta3;
%% for sequential optimization 
Mu=2; lambda = 1; ADTPara = 1; ta =1;
GamIter = zeros(MaxIter,length(regnum)+1); 
Gam_mu_eps = zeros(length(Mu)+length(Epsi)-1,regnum+3);
TempGam = [2 1]; 
for epsi_index = 1:length(Epsi)
iter = 0;
epsi = Epsi(epsi_index)*ones(1,regnum);        
while (abs(TempGam(1)-TempGam(2))>=1e-8)% TempGam(2)<=0.37163   
    iter = iter+1;
    if iter>MaxIter
        break;
    end  
    [Gam,gam,X,Y,Ah,Bh,Ch,Dh,SolverInfo] = multiple_pd_robust(Gasym,TD,d_Thetah,REGID,epsi,XY_PD,Fcn_theta,d_Fcn_theta,FthetaNum,LMI0);       
    iter = iter+1; 
    if iter>MaxIter      
        break;
    end  
end
    mu_epsi_index = mu_index+ epsi_index -1;
    Gam_mu_eps(mu_epsi_index,1) = Gam;
    Gam_mu_eps(mu_epsi_index,2) = mu;
    Gam_mu_eps(mu_epsi_index,3:length(epsi)+2) = epsi;  
    %% check whether all LMIs are satisfied with the results.  
    switch XY_PD
        case 1
            X_theta = X;
            Y_temp = Y;
        case 2
            Y_theta = Y;
            X_temp = X;
        case 3
            X_theta = X;
            Y_theta = Y;
        case 0
            X_temp = X;
            Y_temp = Y;
    end
    lmifail = 0;
end  
gam
[Gam_opt, index] = min(Gam_mu_eps(:,1))
%% recover the controller for time-domain simu.
switch XY_PD 
    case 1       
        X = X_theta;
        Y = Y_temp;
    case 2
        Y = Y_theta;
        X = X_temp;
    case 3
        X = X_theta;
        Y = Y_theta;
    case 0
        X = X_temp;
        Y = Y_temp;
end
