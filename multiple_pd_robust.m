function [Gam,gam,X,Y,Ah,Bh,Ch,Dh,SolverInfo] = multiple_pd_robust(Gasym,TD,d_Thetah,REGID,epsi,XY_PD,Fcn_theta,d_Fcn_theta,FthetaNum, LMI0)
%% This code is for desgin of multiple parameter-dependent robust controller for the miniaturized optical image stabilizers
% Based on the method for switching linear-parameter-varying controller design under uncertain scheduling parameters 
% Modified on 10/02/2015
% created on 10/02/2015

% Gasym, generalized plant, created using Matlab symbols
% ThetaT, cell array, gridded points of theta1 for all regions, ThetaT{1}  for subregion 1, Theta{2} for subregion 2,...
% DeltaT, cell arry, gridded points of delta1 for all regions
% Theta2T, cell arry,gridded points of theta2 for all regions
% Delta2T, cell arry,gridded points of delta2 for all regions
% d_Thetah, cell array, bounds for derivatives of theta, d_Thetah{1} for theta1, d_Thetah{2} for theta2
% REGID, [regid, regid1, regid2]
% epsi, the line search parameter
% XY_PD, 0 for constant X and Y, 1 for 
% Fcn_theta: a function handle, Ftheta(theta) will give all the scalar functions for
% the matrix variables, for instance in X = X0+ f1(theta)X1+f2(theta)X2,
% Fcn_theta(theta) = [1 f1(theta) f2(theta)]
% d_Fcn_theta: a function handle for the derivative of the scalar functions
% FthetaNum: a vector, each element show the number of constant matrices and matrices as a function of GS parameters theta1, theta2, ... 
%            for instance, [1 1 1] means one constant matrix, one matrix as a function of theta1, one matrix as a function of theta2,
%            while [1 2] means one constant matrix, two matrices as a function of theta1, no theta2. 
% LMI0: Initial value for optimization variables when using LMI Lab
%% Paras 
ThetaT = TD.ThetaT;
DeltaT = TD.DeltaT;
Theta2T = TD.Theta2T;
Delta2T = TD.Delta2T;
Theta3T = TD.Theta3T;
Delta3T = TD.Delta3T;
Theta1 = TD.Theta1;
Theta2 = TD.Theta2;
Theta3 = TD.Theta3;
% for avoiding numerical problem
vareps = 1e-5^2; % used for expressing equality constraint
vareps1 = 0; % used for positivity of Lyapunov function
deltam = max(DeltaT{1});
%% 
if sum(abs(Theta3T{1})) == 0 %% only one GS para theta1
   clear Delta3T  % in case 
   Delta3T = {0}; 
   GSParaNum = 2; 
   d_Thetah{3} = 0;
   if max(REGID(:,3)) > 1
       disp('Error!REGID is not for one GS parameter!');
       return;
   end   
else
    GSParaNum = 3;
end
if XY_PD == 0 || XY_PD == 3
    d_Thetah = {0,0,0};
    disp('Derivative of Thetah is set to 0 due to use of constant Lyapunov function or time-invariant scheduling parameter');
end
%% Generalized plant parameter
B2 = Gasym.B2; C2 = Gasym.C2;       
n =  size(Gasym.A,1);
nw = size(Gasym.B1,2);
nu = size(Gasym.B2,2);
nz = size(Gasym.C1,1);
ny = size(Gasym.C2,1); 
I = eye(n); 
%% subregion parameter 
regnum1 = size(ThetaT,2); %num. of subsets for theta1, only partition theta_1
regnum2 = size(Theta2T,2); % num. of subsets for theta2;
regnum = regnum1 * regnum2; % Total num. of subsets. S shape from the bottom to order them, an example:
% 4 5 6
% 1 2 3
% just for positivity of Lyapunov function 
% if GSParaNum == 1
%     theta2gridnum = 1;
% else
%     if length(FthetaNum)<3
%         theta2gridnum = 2;
%     else
%         theta2gridnum = 10;
%     end
% end        
%% Optimization problem definition
lminum = 0;
setlmis([]);
Gam = lmivar(1,[1 1]);
for regid = 1:regnum
    for Id_Ftheta=1:sum(FthetaNum)%
        switch XY_PD
            case 1
                X(Id_Ftheta,regid)=lmivar(1,[n 1]); 
            case 2
                Y(Id_Ftheta,regid)=lmivar(1,[n 1]);
            case 3
                X(Id_Ftheta,regid)=lmivar(1,[n 1]); 
                Y(Id_Ftheta,regid)=lmivar(1,[n 1]);
        end
        Ah(Id_Ftheta,regid)=lmivar(2,[n n]);
        Bh(Id_Ftheta,regid)=lmivar(2,[n ny]);
        Ch(Id_Ftheta,regid)=lmivar(2,[nu n]);
        Dh(Id_Ftheta,regid)=lmivar(2,[nu ny]);
    end
    switch XY_PD
        case 1
            Y(regid)=lmivar(1,[n 1]); 
        case 2
            X(regid)=lmivar(1,[n 1]);
        case 0
           X(regid)=lmivar(1,[n 1]);
           Y(regid)=lmivar(1,[n 1]);
    end            
    gam(regid) = lmivar(1,[1 1]);
end   

%       Ah1 = lmivar(2,[n n]);
%% LMIs for each subregion
for regid = 1:regnum           
   %% Get the grid points for admissible region   
   regid1 = REGID(regid,2); regid2 = REGID(regid,3); 
   thetaT = ThetaT{regid1};deltaT = DeltaT{regid1}; theta2T = Theta2T{regid2}; delta2T = Delta2T{regid2};   theta3T = Theta3T{1};delta3T = Delta3T{1};           
   thetahT = thetaT+deltaT; theta2hT = theta2T + delta2T; theta3hT = theta3T +delta3T;       
%            if UncerType == 1
%                thetahT = thetaT+deltaT; theta2hT = theta2T + delta2T;
%            else
%                thetahT = thetaT.*(1+deltaT); theta2hT = theta2T.*(1 + delta2T);
%            end
   %% Positivity of Lyapunov function
    if XY_PD == 0 % if using PDLF, then there is only one LMI for positivity of LF.
        lminum = lminum + 1;                
        lmiterm([lminum 1 1 X(regid)],1,-1);
        lmiterm([lminum 1 1 0],vareps1); % for making positive definite instead of semi-definite
        lmiterm([lminum 2 1 0],-1);
        lmiterm([lminum 2 2 Y(regid)],1,-1);
        lmiterm([lminum 2 2 0],vareps1);
    else  
        for theta1h  = unique(thetahT)% [min(thetaT),max(thetaT)],  linspace(min(thetaT),max(thetaT),5) % Gridding may not be necessary because use of affine Laypunov function
            for theta2h = unique(theta2hT)%linspace(min(theta2T),max(theta2T),theta2gridnum)     
                for theta3h = unique(theta3hT)
                thetah = [theta1h;theta2h;theta3h];
                Fthetah = Fcn_theta(thetah);                    
                lminum = lminum +1;
                 switch XY_PD 
                     case 1 % PD X                            
                         for Id_Ftheta = 1:length(Fthetah)
                             lmiterm([lminum 1 1 X(Id_Ftheta,regid)],Fthetah(Id_Ftheta),-1);
%                                  X = X+X_theta(:,:,Id_Ftheta,regid)*Fthetah(Id_Ftheta); 
                         end
                         lmiterm([lminum 2 2 Y(regid)],1,-1);
                     case 2 % PD Y
                         for Id_Ftheta = 1:length(Fthetah)
                             lmiterm([lminum 2 2 Y(Id_Ftheta,regid)],Fthetah(Id_Ftheta),-1);
%                                  Y = Y+Y_theta(:,:,Id_Ftheta,regid)*Fthetah(Id_Ftheta); 
                         end   
                         lmiterm([lminum 1 1 X(regid)],1,-1);
                     case 3
                          for Id_Ftheta = 1:length(Fthetah)
                             lmiterm([lminum 1 1 X(Id_Ftheta,regid)],Fthetah(Id_Ftheta),-1);
                             lmiterm([lminum 2 2 Y(Id_Ftheta,regid)],Fthetah(Id_Ftheta),-1);
%                                  X = X+X_theta(:,:,Id_Ftheta,regid)*Fthetah(Id_Ftheta); 
                          end  
                 end                         
                 lmiterm([lminum 1 1 0],vareps1); % for making positive definite instead of semi-definite
                 lmiterm([lminum 2 1 0],-1);    
                 lmiterm([lminum 2 2 0],vareps1);
%                          [X I; I Y] >= eye(2*n)*1e-6; %% This is fairly important for successful simulation in Simulink. If using 0, simulation will be divergent             
                end
            end
        end
    end  

    for Id_theta1 = 1:length(thetaT)
        theta1 = thetaT(Id_theta1);
        for Id_theta2 = 1:length(theta2T)
            theta2 = theta2T(Id_theta2);  
            for Id_theta3 = 1:length(theta3T)
                theta3 = theta3T(Id_theta3); 
                delta1 = deltaT(Id_theta1); delta2 = delta2T(Id_theta2); delta3 = delta3T(Id_theta3);              
%                     if UncerType == 1
%                         theta1h = theta1+delta1; theta2h = theta2+delta2;   
%                     else
%                         theta1h = theta1*(1+delta1); theta2h = theta2*(1+delta2);
%                     end
                theta1h = theta1+delta1; theta2h = theta2+delta2; theta3h = theta3+delta3;
                theta = [theta1;theta2;theta3]; thetah = [theta1h;theta2h;theta3h];        
                Ga = AugPltEval(Gasym, theta);    
                A = Ga.A;
                B1 = Ga.B1;B2 = Ga.B2;
                C1 = Ga.C1;C2 = Ga.C2;
                D11 = Ga.D11; D12 = Ga.D12;
                D21 = Ga.D21; D22 = Ga.D22;       

                Ga =  AugPltEval(Gasym, thetah);
                A_h = Ga.A;  
                Fthetah = Fcn_theta(thetah);   
                d_Fthetah = d_Fcn_theta(thetah);             
                for Id_d_thetah1 = 1:length(d_Thetah{1})
                    d_thetah1 = d_Thetah{1}(Id_d_thetah1);
                    for Id_d_thetah2 = 1:length(d_Thetah{2}) % d_Thetah{2} = 0, for one GS para case
                        d_thetah2 = d_Thetah{2}(Id_d_thetah2); 
                        for Id_d_thetah3 = 1:length(d_Thetah{3})
                            lminum = lminum+1;                            
                            switch XY_PD
                                case 0
                                    lmiterm([lminum 1 1 X(regid)],A,1,'s'); %XA+A'X
                                    lmiterm([lminum 4 1 X(regid)],C1,1);
                                    lmiterm([lminum 2 2 Y(regid)],1,A,'s');
                                    lmiterm([lminum 3 2 Y(regid)],B1',1);                                         
                                case 1
                                    for Id_Ftheta=1:sum(FthetaNum)
                                        lmiterm([lminum 1 1 X(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*A,1,'s');
                                        lmiterm([lminum 4 1 X(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*C1,1);                                        
                                    end
                                    lmiterm([lminum 2 2 Y(regid)],1,A,'s');
                                    lmiterm([lminum 3 2 Y(regid)],B1',1);                                        
                                    % - d_X
                                    for Id_Ftheta = 2:1+FthetaNum(2)
                                        lmiterm([lminum 1 1 X(Id_Ftheta,regid)],-d_Fthetah(Id_Ftheta)*d_thetah1,1);
%                                         -d_X = -(d_X+X_theta(:,:,Id_Ftheta,regid)*d_Fthetah(Id_Ftheta)*d_thetah1); 
                                    end                                 
                                    for Id_Ftheta = 2+FthetaNum(2):sum(FthetaNum)%
                                        lmiterm([lminum 1 1 X(Id_Ftheta,regid)],-d_Fthetah(Id_Ftheta)*d_thetah2,1);
                                    end
                                case 2
                                    for Id_Ftheta = 1:sum(FthetaNum)
                                        lmiterm([lminum 2 2 Y(Id_Ftheta,regid)],Fthetah(Id_Ftheta),A,'s');
                                        lmiterm([lminum 3 2 Y(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*B1',1);  
                                    end
                                    lmiterm([lminum 1 1 X(regid)],A,1,'s'); %XA+A'X
                                    lmiterm([lminum 4 1 X(regid)],C1,1);   
                                    % d_Y
                                    for Id_Ftheta = 2:1+FthetaNum(2)
                                        lmiterm([lminum 2 2 Y(Id_Ftheta,regid)],d_Fthetah(Id_Ftheta)*d_thetah1,1);
%                                         d_X = d_X+X_theta(:,:,Id_Ftheta,regid)*d_Fthetah(Id_Ftheta)*d_thetah1; 
                                    end                                 
                                    for Id_Ftheta = 2+FthetaNum(2):sum(FthetaNum)% if 
                                        lmiterm([lminum 2 2 Y(Id_Ftheta,regid)],d_Fthetah(Id_Ftheta)*d_thetah2,1);
                                    end
                                case 3
                                     for Id_Ftheta = 1:sum(FthetaNum)
                                        lmiterm([lminum 1 1 X(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*A,1,'s');
                                        lmiterm([lminum 4 1 X(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*C1,1);   
                                        lmiterm([lminum 2 2 Y(Id_Ftheta,regid)],Fthetah(Id_Ftheta),A,'s');
                                        lmiterm([lminum 3 2 Y(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*B1',1);  
                                     end                                   
                            end   
                            lmiterm([lminum 2 1 0],A');
                            lmiterm([lminum 3 1 0],B1');
                            lmiterm([lminum 3 3 gam(regid)],-1,1);
                            lmiterm([lminum 4 2 0],C1);
                            lmiterm([lminum 4 3 0],D11);
                            lmiterm([lminum 4 4 gam(regid)],-1,1);   

                            for Id_Ftheta=1:sum(FthetaNum)   
%                                 regidtemp = regid;
%                                 regid = 2; 
                                lmiterm([lminum 2 1 Ah(Id_Ftheta,regid)],Fthetah(Id_Ftheta),1);                                  
                                lmiterm([lminum 2 2 Bh(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*1,C2,'s'); 
                                lmiterm([lminum 3 2 -Bh(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*D21',1);
                                lmiterm([lminum 1 1 Ch(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*B2,1,'s');
                                lmiterm([lminum 4 1 Ch(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*D12,1);                                
                                lmiterm([lminum 2 1 -Dh(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*C2',B2');                                
                                lmiterm([lminum 3 1 -Dh(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*D21',B2'); 
                                lmiterm([lminum 4 2 Dh(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*D12,C2); 
                                lmiterm([lminum 4 3 Dh(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*D12,D21);    
%                                 regid = regidtemp;
                            end   
                            % for uncertainty
                            if deltam ~= 0
                                switch XY_PD
                                    case 0                                        
                                        lmiterm([lminum 5 1 X(regid)],epsi(regid),1);
                                        lmiterm([lminum 5 2 Y(regid)],(A-A_h)',1);
                                        lmiterm([lminum 5 5 0],-epsi(regid));                                            
                                    case 1                                        
                                        for Id_Ftheta=1:sum(FthetaNum)
                                            lmiterm([lminum 5 1 X(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*epsi(regid),1); %XA+A'X
                                        end
                                        lmiterm([lminum 5 2 Y(regid) ],(A-A_h)',1);
                                        lmiterm([lminum 5 5 0],-epsi(regid));
                                    case 2
                                        lmiterm([lminum 5 1 X(regid)],epsi(regid),1);
                                        for Id_Ftheta=1:sum(FthetaNum)
                                            lmiterm([lminum 5 2 Y(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*(A-A_h)',1); %XA+A'X
                                        end
                                        lmiterm([lminum 5 5 0],-epsi(regid)); 
                                    case 3
                                        for Id_Ftheta=1:sum(FthetaNum)
                                            lmiterm([lminum 5 1 X(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*epsi(regid),1);
                                            lmiterm([lminum 5 2 Y(Id_Ftheta,regid)],Fthetah(Id_Ftheta)*(A-A_h)',1); %XA+A'X
                                        end
                                        lmiterm([lminum 5 5 0],-epsi(regid));                                               
                                end 
                            end %deltam 
                        end %Id_d_thetah3
                    end %Id_d_thetah2
                end %Id_d_thetah1
            end %Id_theta3
        end % Id_theta2
    end % Id_theta1
    lminum = lminum + 1;
    lmiterm([lminum 1 1 gam(regid)],1,1);
    lmiterm([lminum 1 1 Gam],-1,1);
end %regid    

lmisys = getlmis;
 nvar = decnbr(lmisys); %number of decision variable.
c = zeros(nvar,1);
c(1)=1; 
options(1)= 1e-4; % relative accuary on the optimal value
options(2)= 400; %Number of Iteration
options(4) =  10; % J, the code terminates when the objective has not decreased by more than the desired relative accuracy during the last J iterations
options(5)= 0; % 1 for not showing the process
% specifying the initial value for the decision variables
% if ~isempty(LMI0)  
%     if regnum == 1 && XY_PD == 1 && sum(FthetaNum) == 3                 
%          xinit = mat2dec(lmisys,LMI0.Gam,...
%             LMI0.X(:,:,1),LMI0.Ah(:,:,1),LMI0.Bh(:,:,1),LMI0.Ch(:,:,1),LMI0.Dh(:,:,1),...
%             LMI0.X(:,:,2),LMI0.Ah(:,:,2),LMI0.Bh(:,:,2),LMI0.Ch(:,:,2),LMI0.Dh(:,:,2),...
%             LMI0.X(:,:,3),LMI0.Ah(:,:,3),LMI0.Bh(:,:,3),LMI0.Ch(:,:,3),LMI0.Dh(:,:,3),...
%             LMI0.Gam, ...
%             LMI0.Y);
%     elseif regnum == 2
%         switch XY_PD
%             case 1
%                 if sum(FthetaNum) == 2
%                     xinit = mat2dec(lmisys,LMI0.Gam,....
%                         LMI0.X(:,:,1),LMI0.Ah(:,:,1),LMI0.Bh(:,:,1),LMI0.Ch(:,:,1),LMI0.Dh(:,:,1),...
%                         LMI0.X(:,:,2),LMI0.Ah(:,:,2),LMI0.Bh(:,:,2),LMI0.Ch(:,:,2),LMI0.Dh(:,:,2),...
%                         LMI0.Gam,...
%                         LMI0.X(:,:,1),LMI0.Ah(:,:,1),LMI0.Bh(:,:,1),LMI0.Ch(:,:,1),LMI0.Dh(:,:,1),...
%                         LMI0.X(:,:,2),LMI0.Ah(:,:,2),LMI0.Bh(:,:,2),LMI0.Ch(:,:,2),LMI0.Dh(:,:,2),...
%                         LMI0.Gam,...
%                         LMI0.Y);
%                 elseif sum(FthetaNum) == 3
%                      xinit = mat2dec(lmisys,LMI0.Gam,...
%                         LMI0.X(:,:,1),LMI0.Ah(:,:,1),LMI0.Bh(:,:,1),LMI0.Ch(:,:,1),LMI0.Dh(:,:,1),...
%                         LMI0.X(:,:,2),LMI0.Ah(:,:,2),LMI0.Bh(:,:,2),LMI0.Ch(:,:,2),LMI0.Dh(:,:,2),...
%                         LMI0.X(:,:,3),LMI0.Ah(:,:,3),LMI0.Bh(:,:,3),LMI0.Ch(:,:,3),LMI0.Dh(:,:,3),...
%                         LMI0.Gam, ...
%                         LMI0.X(:,:,1),LMI0.Ah(:,:,1),LMI0.Bh(:,:,1),LMI0.Ch(:,:,1),LMI0.Dh(:,:,1),...
%                         LMI0.X(:,:,2),LMI0.Ah(:,:,2),LMI0.Bh(:,:,2),LMI0.Ch(:,:,2),LMI0.Dh(:,:,2),...
%                         LMI0.X(:,:,3),LMI0.Ah(:,:,3),LMI0.Bh(:,:,3),LMI0.Ch(:,:,3),LMI0.Dh(:,:,3),...
%                         LMI0.Gam,...
%                         LMI0.Y);
%                 end                    
%             case 2
%                 if sum(FthetaNum) == 2
%                     xinit = mat2dec(lmisys,LMI0.Gam,...
%                         LMI0.Y(:,:,1),LMI0.Ah(:,:,1),LMI0.Bh(:,:,1),LMI0.Ch(:,:,1),LMI0.Dh(:,:,1),...
%                         LMI0.Y(:,:,2),LMI0.Ah(:,:,2),LMI0.Bh(:,:,2),LMI0.Ch(:,:,2),LMI0.Dh(:,:,2),...
%                         LMI0.Gam,...
%                         LMI0.Y(:,:,1),LMI0.Ah(:,:,1),LMI0.Bh(:,:,1),LMI0.Ch(:,:,1),LMI0.Dh(:,:,1),...
%                         LMI0.Y(:,:,2),LMI0.Ah(:,:,2),LMI0.Bh(:,:,2),LMI0.Ch(:,:,2),LMI0.Dh(:,:,2),...
%                         LMI0.Gam,...
%                         LMI0.X);                        
%                 elseif sum(FthetaNum) == 3
%                     xinit = mat2dec(lmisys,LMI0.Gam,...
%                         LMI0.Y(:,:,1),LMI0.Ah(:,:,1),LMI0.Bh(:,:,1),LMI0.Ch(:,:,1),LMI0.Dh(:,:,1),...
%                         LMI0.Y(:,:,2),LMI0.Ah(:,:,2),LMI0.Bh(:,:,2),LMI0.Ch(:,:,2),LMI0.Dh(:,:,2),...
%                         LMI0.Y(:,:,3),LMI0.Ah(:,:,3),LMI0.Bh(:,:,3),LMI0.Ch(:,:,3),LMI0.Dh(:,:,3),...
%                         LMI0.Gam, ...
%                         LMI0.Y(:,:,1),LMI0.Ah(:,:,1),LMI0.Bh(:,:,1),LMI0.Ch(:,:,1),LMI0.Dh(:,:,1),...
%                         LMI0.Y(:,:,2),LMI0.Ah(:,:,2),LMI0.Bh(:,:,2),LMI0.Ch(:,:,2),LMI0.Dh(:,:,2),...
%                         LMI0.Y(:,:,3),LMI0.Ah(:,:,3),LMI0.Bh(:,:,3),LMI0.Ch(:,:,3),LMI0.Dh(:,:,3),...
%                         LMI0.Gam,...
%                         LMI0.X);
%                 end
%             case 0
%                 if sum(FthetaNum) == 2
%                     xinit = mat2dec(lmisys,LMI0.Gam,...
%                         LMI0.Ah(:,:,1),LMI0.Bh(:,:,1),LMI0.Ch(:,:,1),LMI0.Dh(:,:,1),...
%                         LMI0.Ah(:,:,2),LMI0.Bh(:,:,2),LMI0.Ch(:,:,2),LMI0.Dh(:,:,2),...
%                         LMI0.Gam,...
%                         LMI0.Ah(:,:,1),LMI0.Bh(:,:,1),LMI0.Ch(:,:,1),LMI0.Dh(:,:,1),...
%                         LMI0.Ah(:,:,2),LMI0.Bh(:,:,2),LMI0.Ch(:,:,2),LMI0.Dh(:,:,2),...
%                         LMI0.Gam,...
%                         LMI0.Y,LMI0.Y);
%                 elseif sum(FthetaNum) == 3
%                      xinit = mat2dec(lmisys,LMI0.Gam,...
%                         LMI0.Ah(:,:,1),LMI0.Bh(:,:,1),LMI0.Ch(:,:,1),LMI0.Dh(:,:,1),...
%                         LMI0.Ah(:,:,2),LMI0.Bh(:,:,2),LMI0.Ch(:,:,2),LMI0.Dh(:,:,2),...
%                         LMI0.Ah(:,:,3),LMI0.Bh(:,:,3),LMI0.Ch(:,:,3),LMI0.Dh(:,:,3),...
%                         LMI0.Gam,...
%                         LMI0.Ah(:,:,1),LMI0.Bh(:,:,1),LMI0.Ch(:,:,1),LMI0.Dh(:,:,1),...
%                         LMI0.Ah(:,:,2),LMI0.Bh(:,:,2),LMI0.Ch(:,:,2),LMI0.Dh(:,:,2),...
%                         LMI0.Ah(:,:,3),LMI0.Bh(:,:,3),LMI0.Ch(:,:,3),LMI0.Dh(:,:,3),...
%                         LMI0.Gam,...
%                         LMI0.X,LMI0.Y);
%                 end
%         end
%     [Gam,xopt] = mincx(lmisys,c,options,xinit);
%     end
% end

[Gam,xopt] = mincx(lmisys,c,options);


%% check whether all constraints are satisfied
%         lmifail = 0;
%         evals = evallmi(lmisys,xopt);        
%         for i = 1: lminum
%             [lhs,rhs] = showlmi(evals,i);
%             if max(real(eig(lhs-rhs))) > 5e-7 
%                 eig_max = max(real(eig(lhs-rhs)))
%                 disp ('Not all LMI constranits are satisfied')
%                 lmifail = 1;
% %                 break;
%             end
%         end 
%         lmifail
if ~isequal(Gam,[])       
    for regid = 1:regnum
    for Id_Ftheta = 1:sum(FthetaNum)
        Akh(:,:,Id_Ftheta,regid) = dec2mat(lmisys,xopt,Ah(Id_Ftheta,regid));
        Bkh(:,:,Id_Ftheta,regid) = dec2mat(lmisys,xopt,Bh(Id_Ftheta,regid));
        Ckh(:,:,Id_Ftheta,regid) = dec2mat(lmisys,xopt,Ch(Id_Ftheta,regid));
        Dkh(:,:,Id_Ftheta,regid) = dec2mat(lmisys,xopt,Dh(Id_Ftheta,regid));
        switch XY_PD
            case 1
                Xk(:,:,Id_Ftheta,regid) = dec2mat(lmisys,xopt,X(Id_Ftheta,regid));
            case 2
                Yk(:,:,Id_Ftheta,regid) = dec2mat(lmisys,xopt,Y(Id_Ftheta,regid)); 
            case 3
                Xk(:,:,Id_Ftheta,regid) = dec2mat(lmisys,xopt,X(Id_Ftheta,regid));
                Yk(:,:,Id_Ftheta,regid) = dec2mat(lmisys,xopt,Y(Id_Ftheta,regid)); 
        end
    end
    switch XY_PD 
        case 1
            Yk(:,:,regid) = dec2mat(lmisys,xopt,Y(regid));
        case 2
            Xk(:,:,regid) = dec2mat(lmisys,xopt,X(regid));
        case 0
            Xk(:,:,regid) = dec2mat(lmisys,xopt,X(regid));
            Yk(:,:,regid) = dec2mat(lmisys,xopt,Y(regid));
    end
    gam_reg(regid)  = dec2mat(lmisys,xopt,gam(regid));
    end 
    Ah = Akh;
    Bh = Bkh;
    Ch = Ckh;
    Dh = Dkh;        
    X = Xk;
    Y = Yk;   
    gam = gam_reg;
    lmifail = 0;
else
    Ah = [];
    Bh = [];
    Ch = [];
    Dh = [];        
    X = [];
    Y = [];   
    gam = [];   
    lmifail = 1
    Gam = NaN;
end      
clear Akh Bkh Ckh Dkh Xk Yk     
SolverInfo.lmifail = lmifail;
