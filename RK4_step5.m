%% RK4 for step 5

function x_hat = RK4_step5(A,B,C,x0,K,L,T,x)
t=0:0.01:T;
   x_hat=zeros(6,length(t));
   x_hat(:,1)=x0;
   h=0.01;
   for i=1:length(t)-1
        f1 = RHS_obs(x_hat(:, i), A, B,C,K,L,x(:,i));
        f2 = RHS_obs(x_hat(:, i)+f1/2*h, A, B,C,K,L,x(:,i));
        f3 = RHS_obs(x_hat(:, i)+f2/2*h, A, B,C,K,L,x(:,i));
        f4 = RHS_obs(x_hat(:, i)+f3*h, A, B,C,K,L,x(:,i));
        x_hat(:, i+1) = x_hat(:, i) + h*(f1/6+((f2+f3)/3)+f4/6);
   end

   function dx_hat=RHS_obs(x_hat,A,B,C,K,L,x)
        dx_hat = (A-L*C)*x_hat - B*K*x_hat + L*C*x;
   end
end