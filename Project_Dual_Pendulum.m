%%
% https://github.com/Karopsy/Dual_Pendulum.git
%% PHASE A - Part 1 : Swing-up
% STEP 1 : Sworm training : (with optimal s.QT matrix to compute the cost function)
close all
clear u_k
clear x_k
load s_struct.mat
target_ang=deg2rad(15);
target_vel=deg2rad(30);
u_k(:,1)=zeros(s.N+1,1);
for i=1:50,i
    [u_k(:,i+1),x_k]=Example_22_2(s.T,u_k(:,i));
    figure(4); 
    grid minor
    scatter(i,rad2deg(x_k(2,end)),'r')
    hold on 
    scatter(i,rad2deg(x_k(3,end)),'b')
    scatter(i,rad2deg(x_k(5,end)),'g')
    scatter(i,rad2deg(x_k(6,end)),'k')
    yline([-rad2deg(target_ang) rad2deg(target_ang)])
    yline([-rad2deg(target_vel) rad2deg(target_vel)])
    legend("\Theta_1","\Theta_2","Ang vel \Theta_1","Ang vel \Theta_2")
    ylabel("Angles - deg")
    xlabel('Iterations')
    if -target_ang<x_k(2,end) && x_k(2,end)<target_ang && -target_ang<x_k(3,end) && x_k(3,end)<target_ang && abs(x_k(5,end))<target_vel && abs(x_k(6,end))<target_vel
        break
    end
end

%% Phase A - Part 2 : Optimal Control LTV (K(t))
Q=eye(6);
QT=eye(6);%diag([5 20 20 0.1 60 30]);%diag([5 20 20 0.1 60 30]);
R=5;

for i =1:length(x_k)
    A(:,:,i)=Compute_A(x_k(:,i),s);
    E(:,:,i)=Compute_E(x_k(:,i),s);
    B(:,:,i)=[0 0 0 1 0 0]';
end    
for i=length(x_k):-1:2 %Solve it backward in time
    if i==length(x_k)
        X(:,:,i)=QT;
    end
    f1=DRE(X(:,:,i),A(:,:,i),B(:,:,i),R,E(:,:,i),Q); 
    f2=DRE(X(:,:,i)+s.h*f1/2,A(:,:,i),B(:,:,i),R,E(:,:,i),Q);    
    f3=DRE(X(:,:,i)+s.h*f2/2,A(:,:,i),B(:,:,i),R,E(:,:,i),Q); 
    f4=DRE(X(:,:,i)+s.h*f3,A(:,:,i),B(:,:,i),R,E(:,:,i),Q);
    X(:,:,i-1)=X(:,:,i)+s.h*(f1/6+(f2+f3)/3+f4/6); 
    K(:,i)=-R^-1*B(:,:,i)'*X(:,:,i)*E(:,:,i);
    EV_contr(:,i)=eig(E(:,:,i)^-1*A(:,:,i)+E(:,:,i)^-1*B(:,:,i)*K(:,i)');
end
%% Verification of K(t)
% Because EV_contr are not neg then we are checking that the x_prim path
% is close to the optimized path:

for i=1:length(x_k)-1
    if i==1
        x_prim(:,1)=x_k([1:6],1);
    end
    f1=RHS_modi(x_prim(:,i),A(:,:,i),B(:,:,i),E(:,:,i),K(:,i)',u_k(i,end));
    f2=RHS_modi(x_prim(:,i)+s.h*f1/2,A(:,:,i),B(:,:,i),E(:,:,i),K(:,i)',u_k(i,end));   
    f3=RHS_modi(x_prim(:,i)+s.h*f2/2,A(:,:,i),B(:,:,i),E(:,:,i),K(:,i)',u_k(i,end)); 
    f4=RHS_modi(x_prim(:,i)+s.h*f3,A(:,:,i),B(:,:,i),E(:,:,i),K(:,i)',u_k(i,end));
    x_prim(:,i+1)=x_prim(:,i)+s.h*(f1/6+(f2+f3)/3+f4/6);
    x_sim(:,i)=x_k([1:6],i)+s.h*x_prim(:,i);
end
x_sim(:,length(x_k))=x_k([1:6],end)+s.h*x_prim(:,end);
%%
figure(10)
sgtitle("System with Optimal Control - LTV ")
subplot(3,2,1)
plot(s.t,x_k(1,:))
hold on 
plot(s.t,x_sim(1,:))
xlabel("Time (sec)")
ylabel("Position of the cart - m")
legend("Guide","Opti Contr Traj")

subplot(3,2,2)
plot(s.t,x_k(2,:))
hold on 
plot(s.t,x_sim(2,:))
xlabel("Time (sec)")
ylabel("Angle \Theta_1 (rad)")
legend("Guide","Opti Contr Traj")

subplot(3,2,3)
plot(s.t,x_k(3,:))
hold on 
plot(s.t,x_sim(3,:))
xlabel("Time (sec)")
ylabel("Angle \Theta_2 (rad)")
legend("Guide","Opti Contr Traj")

subplot(3,2,4)
plot(s.t,x_k(4,:))
hold on 
plot(s.t,x_sim(4,:))
xlabel("Time (sec)")
ylabel("Velocity of the cart m/s")
legend("Guide","Opti Contr Traj")

subplot(3,2,5)
plot(s.t,x_k(5,:))
hold on 
plot(s.t,x_sim(5,:))
xlabel("Time (sec)")
ylabel("Angular velocity \Theta_1 (rad/s)")
legend("Guide","Opti Contr Traj")

subplot(3,2,6)
plot(s.t,x_k(6,:))
hold on 
plot(s.t,x_sim(6,:))
xlabel("Time (sec)")
ylabel("Angular velocity \Theta_2 (rad/s)")
legend("Guide","Opti Contr Traj")



%% Phase A - Part 3 : Kalman Filter LTV (L(t))
Q=eye(6);
P0=eye(6);%diag([5 20 20 0.1 60 30]);
R=diag([5 5 5 5 5 5]);

for i=1:length(x_k)
    C(:,:,i)=diag([1 1 1 0 0 0]);
    E_kal(:,:,i)=eye(6,6);
end

for i=1:length(x_k)-1 %Solve it forward in time
    if i==1
        P(:,:,i)=P0;
    end
    f1=DRE(P(:,:,i),A(:,:,i)',C(:,:,i)',R,E_kal(:,:,i),Q); 
    f2=DRE(P(:,:,i)+s.h*f1/2,A(:,:,i)',C(:,:,i)',R,E_kal(:,:,i),Q);   
    f3=DRE(P(:,:,i)+s.h*f2/2,A(:,:,i)',C(:,:,i)',R,E_kal(:,:,i),Q); 
    f4=DRE(P(:,:,i)+s.h*f3,A(:,:,i)',C(:,:,i)',R,E_kal(:,:,i),Q);
    P(:,:,i+1)=P(:,:,i)+s.h*(f1/6+(f2+f3)/3+f4/6); 
    %K(:,i)=-R^-1*B'*X(:,:,i)*E(:,:,i);
    L(:,:,i)=-P(:,:,i)*C(:,:,i)'*R^-1;
    EV_obs(:,:,i)=eig(A(:,:,i)+L(:,:,i)*C(:,:,i));

end

L(:,:,length(x_k))=-P(:,:,length(x_k))*C(:,:,length(x_k)-1)'*R^-1; %C(:,:,length(x_k)-1) = C(:,:,length(x_k)) C is constant
EV_obs(:,:,length(x_k))=eig(A(:,:,length(x_k))+L(:,:,length(x_k))*C(:,:,length(x_k)-1));
%% Verification of L(t)
% Because EV_obs are not neg then we do the same thing as for the
% controller : checking that the x_hat path is close to the optimized path:

x_k=x_k([1:6],:);
for i=1:length(x_k)-1
    if i==1
        %x_hat(:,1)=x_k([1:6],1);
        x_hat(:,1)=zeros(6,1);
    end
    f1=RHS_kalman(x_hat(:,i),E(:,:,i)^-1*A(:,:,i),L(:,:,i),C(:,:,i),E(:,:,i)^-1*B(:,:,i),K(:,i)',x_k(:,end),x_prim(:,i),u_k(i,end));
    f2=RHS_kalman(x_hat(:,i)+s.h*f1/2,E(:,:,i)^-1*A(:,:,i),L(:,:,i),C(:,:,i),E(:,:,i)^-1*B(:,:,i),K(:,i)',x_k(:,end),x_prim(:,i),u_k(i,end));   
    f3=RHS_kalman(x_hat(:,i)+s.h*f2/2,E(:,:,i)^-1*A(:,:,i),L(:,:,i),C(:,:,i),E(:,:,i)^-1*B(:,:,i),K(:,i)',x_k(:,end),x_prim(:,i),u_k(i,end)); 
    f4=RHS_kalman(x_hat(:,i)+s.h*f3,E(:,:,i)^-1*A(:,:,i),L(:,:,i),C(:,:,i),E(:,:,i)^-1*B(:,:,i),K(:,i)',x_k(:,end),x_prim(:,i),u_k(i,end));
    x_hat(:,i+1)=x_hat(:,i)+s.h*(f1/6+(f2+f3)/3+f4/6);
    x_sim2(:,i)=x_k(:,i)+s.h*x_hat(:,i);
end
x_sim2(:,length(x_k))=x_k(:,end)+s.h*x_hat(:,end);
%% Figures

figure(20)
sgtitle("System with Kalman Filter - LTV ")
subplot(3,2,1)
plot(s.t,x_k(1,:))
hold on 
plot(s.t,x_sim2(1,:))
xlabel("Time (sec)")
ylabel("Position of the cart - m")
legend("Guide","Kalman Traj")

subplot(3,2,2)
plot(s.t,x_k(2,:))
hold on 
plot(s.t,x_sim2(2,:))
xlabel("Time (sec)")
ylabel("Angle \Theta_1 (rad)")
legend("Guide","Kalman Traj")

subplot(3,2,3)
plot(s.t,x_k(3,:))
hold on 
plot(s.t,x_sim2(3,:))
xlabel("Time (sec)")
ylabel("Angle \Theta_2 (rad)")
legend("Guide","Kalman Traj")

subplot(3,2,4)
plot(s.t,x_k(4,:))
hold on 
plot(s.t,x_sim2(4,:))
xlabel("Time (sec)")
ylabel("Velocity of the cart m/s")
legend("Guide","Kalman Traj")

subplot(3,2,5)
plot(s.t,x_k(5,:))
hold on 
plot(s.t,x_sim2(5,:))
xlabel("Time (sec)")
ylabel("Angular velocity \Theta_1 (rad/s)")
legend("Guide","Kalman Traj")

subplot(3,2,6)
plot(s.t,x_k(6,:))
hold on 
plot(s.t,x_sim2(6,:))
xlabel("Time (sec)")
ylabel("Angular velocity \Theta_2 (rad/s)")
legend("Guide","Kalman Traj")



%% Phase B : Steps 4 and 5 (Stabilization)

% Creating state space system to solve for K and L
g=9.8;
E=Compute_E(x_k(1:6,end),s);
N=[0,0,0,1,0,0;
    0,0,0,0,1,0;
    0,0,0,0,0,1;
    0,0,0,0,0,0,;
    0,s.m1*g*s.L1,0,0,0,0;
    0,0,s.m2*g*s.L2,0,0,0];

A=E^-1*N;
B=E^-1*[0,0,0,1,0,0]';
C=diag([1 1 1 0 0 0]);

D = zeros(size(C,1),size(B,2));

Qc=diag([100 100 1 1 100 100]);
Qe=eye(6);
Rc=3;
Re=9;
sys=ss(A,B,C,D);
[X_b,K_b,EV_contr_b]=icare(A,B,Qc,Rc); 
[P_b,L_b,EV_obs_b]=icare(A',C,Qe,Re);

sysKF = ss(A-L_b'*C,[B L_b'],eye(6),0*[B L_b']);  % Kalman filter estimator

T2=15;
t=0:0.01:T2;

x=zeros(6,length(t));
x0=x_k(1:6,end);
%x0=[0.1 0.1 0.1 0.1 0.1 0.1];
x = RK4_step4(A,B,x0,K_b,T2);

figure(21)
for  i=1:6
    hold on
    plot(t,x(i,:))
end
xlabel("Time (sec)")
ylabel("State Variables")
legend("Cart pos","Angle \Theta_1",'Angle \Theta_2',"Cart vel","Ang vel 1","Ang vel 2")
title("System at \infty horizon with Optimal Control ")

x_hat=zeros(6,length(t));
x_hat = RK4_step5(A,B,C,x0,K_b,L_b',T2,x);

figure(31)
for i=1:6
    hold on
    plot(t,flip(x_hat(i,:)))
end
xlabel("Time (sec)")
ylabel("State Variables")
legend("Cart pos","Angle \Theta_1",'Angle \Theta_2',"Cart vel","Ang vel 1","Ang vel 2")
title("System at \infty horizon with Kalman Filter ")

