
%% obtain the LTI controllers from the parameter-dependent robust controllers
% and anlayse the stability of the closed-loop system
% close all;

% simin =[0 0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008; -1 -0.8 -0.5 -0.2 0 0.2 0.5 0.8 1]';
% simin =[0 1e-3:1e-3:4e-3; linspace(Theta1(1), Theta1(end),5)]';
W2n = [97.11 101.6 101.7 104.9 106.8]*2*pi;
W1n = [61.556 62.84 63.84 67.33 68.675]*2*pi;
% simin = [0 1:length(W2n)-1; W2n]';
% simtime = size(simin,1)-1;
B2 = Gasym.B2; C2 = Gasym.C2;
if SwLogic == 0
    REGid = ones(1,5);
else 
    REGid = [1 1 1 4 4];
end

%% Get the five controllers
for i = 1:length(W2n)
    theta1 = W2n(i);
    theta2 = W1n(i);
    theta3 = 0; 
    A_t = [                   0,       theta2,                 0,             0,       0,     0,       0;
             -1.0*theta2, -0.01*theta2,                 0,             0,    54.4,     0,       0;
                       0,            0,                 0,         639.0,       0,     0,       0;
                       0,            0, -0.00156*theta1^2, -0.001*theta1,  2188.0,     0,       0;
                       0,            0,                 0,             0, -1744.0,     0,       0;
     -21.3*theta2*theta3,            0,           -1400.0,             0,       0, -19.8,       0;
     -64.0*theta2*theta3,            0,           -4211.0,             0,       0,     0, -5000.0];
    reg = [Theta1;0 0;0 0];
    CtrlSel = REGid(i);
%     if SwLogic == 0
%         CtrlSel = 1;    
%     else
%         if theta1 > reg(1,1) && theta1 <= reg(1,2)
%            CtrlSel = 1;
%         elseif  theta1 > reg(2,1) && theta1 <= reg(2,2) %subset 2
%             CtrlSel = 2;
%         elseif theta1 >  reg(3,1) && theta1 <= reg(3,2) 
%             CtrlSel = 3;
%         end
%     end   
    regid = CtrlSel;
    Ah_t = Ah(:,:,1,regid)+theta1*Ah(:,:,2,regid)+theta2*Ah(:,:,3,regid);
    Bh_t = Bh(:,:,1,regid)+theta1*Bh(:,:,2,regid)+theta2*Bh(:,:,3,regid);
    Ch_t = Ch(:,:,1,regid)+theta1*Ch(:,:,2,regid)+theta2*Ch(:,:,3,regid);
    Dh_t = Dh(:,:,1,regid)+theta1*Dh(:,:,2,regid)+theta2*Dh(:,:,3,regid);
    regid)+theta(1)*Dh(:,:,2,regid)+theta(1)^2*Dh(:,:,3,regid);

    if XY_PD == 3
        Y_t = Y(:,:,1,regid)+theta1*Y(:,:,2,regid)+theta2*Y(:,:,3,regid);
        X_t = X(:,:,1,regid)+theta1*X(:,:,2,regid)+theta2*X(:,:,3,regid);
        M = -X_t;
        N =  Y_t - eye(n)/X_t; 
    %     Minv = M\eye(n);
    %     Ninv = (Y_t + Minv)\eye(n);
    else
        X_t = X(:,:,1);
        Y_t = Y(:,:,1);
        M = -X_t;
        N =  Y_t - eye(n)/X_t;   
    end
    clear Dki Cki Bki Aki
    Dk = Dh_t ;
    Ck = (Ch_t - Dk*C2*X_t)/M';
    Bk = N\(Bh_t -Y_t*B2*Dk);
    Ak =  N\(Ah_t -Bh_t*C2*X_t -...
        Y_t*B2*Ch_t -Y_t*(A_t-B2*Dh_t*C2)*X_t)/M';  
        
    Kmgs{i} = ss(Ak,Bk,Ck,Dk)*Fy;
end
% close all; 
%% Model reduction to improve numerical stability
% opt = balredOptions('StateElimMethod','Truncate');
% figure;
% for i = 1: 5
%     Kgsr{i} = balred(Kgs{i},order(Kgs{i})-2); bode(Kgs{i},Kgsr{i}); cond(Kgsr{i}.a)
% end
%% Discretize the controller and expressed in transfer function
figure;
Ts = 5e-4;
opt = c2dOptions('Method','tustin','PrewarpFrequency',70*2*pi);
for i = 1: 5   
    Kpdrd{i} = c2d(tf(Kpdr{i}),Ts,opt); bode(Kpdr{i},Kpdrd{i});hold on;
    [Kpdrd_num{i},Kpdrd_den{i}] = tfdata(Kpdrd{i},'v'); 
end
return;
%% Paras for plant
% DCgain = 10^(12.29/20); 
% zeta1 = 5e-3; 
% w10 = 63.84*2*pi;
% zeta2 = 5e-4;
% Td = 21.86/(360*52.87); % pure time delay
% % Theta1 = [96 108]*2*pi;
%%
hf1 = figure('Name','Bode plot of plant and controller'); 
for k=1:2
    h(k) = subplot(2,1,k);
end
hf2 = figure('Name','Bode plot of open loop system');
hf3 = figure('Name','Sensitivity , 1/We, 1/Wu');
h3 = subplot(2,1,1);
h4 = subplot(2,1,2);
% hf4 = figure('Name','Bode plot of controller');
% hf5 = figure('Name','Complementary sensitivity');
deltam = 0;

K1 = 0.18;
K2 = DCgain;
for i = 1:size(W2n,1)
    theta1 = W2n(i);
    
    A1 = [0 w10/w20*theta1; -w10/w20*theta1 -2*zeta1*w10/w20*theta1]; %use 1 instead of /w20*theta1 to make B term constant
    B1 = [0 1]';
    C1 = [K1*w10*1/w20*theta1 0];
    D1 = 0;

    %% 2nd mode
    A2 = [0 w20;-theta1^2/w20 -2*zeta2*theta1];% use w20 instead of theta1 to make B term constant
    B2 = [0 1]'*40;
    C2 = [K2*w20 0]/40;
    D2 = 0;
 
    % pade approximation of pure delay
    tf3 = ss(tf([-0.5*Td  1],[0.5*Td  1]));
    A3 = tf3.a; B3 = tf3.b; C3 = tf3.c; D3 = tf3.d;

    % sum of (A1,B1,C1,D1) and (A2,B2,C2,D2)
    A = [A1 zeros(size(A1,1),size(A2,2)); zeros(size(A2,1),size(A1,2)) A2];
    B = [B1;B2];
    C = [C1 C2];
    D = D1+D2;

    % product of (A,B,C,D) and (A3,B3,C3,D3)
    A = [A B*C3; zeros(size(A3,1),size(A,2)) A3];
    B = [B*D3; B3];
    C = [C D*C3];
    D = D*D3;
        
    Plt = ss(A,B,C,D);    
    
    Ctrk =  Kpdrr{1}; 
%     Ctrk = sysKmuR;
    OpenLoop = Ctrk*Plt;
    figure(hf1);
    subplot(h(1))
    bodemag(Plt,{1e-1,1e5});hold on;
    subplot(h(2));
%     figure(hf4); 
    bodemag(Ctrk,{1e-1,1e5});hold on;   
    figure(hf2);
    bodemag(OpenLoop,{1e-1,1e5});hold on;     
    %% close loop stability analysis
    sys = feedback(1,OpenLoop);    
    figure(hf3);
    subplot(h3);
    bodemag(sys,{1e-1,1e5});hold on;
    subplot(h4);
    bodemag(sys*Ctrk,{1e-1,1e5});hold on;
    eigs = eig(sys.a);
    if max(real(eigs)) >= 0
        disp('Closed loop system is unstable at');
        theta1
        deltam
        Eig = eigs
    end
     kk = 1;
     %% complementary sensitivity
%      figure(hf5);
%      bodemag(1-sys,{1e-1,1e3});hold on;
end
figure(hf3);
subplot(h3);
bodemag(1/We);
subplot(h4);
bodemag(ss(1/Wu));
%%