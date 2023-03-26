%% RK4 for step 4

function x = RK4_step4(A,B,x0,K,T)
t=0:0.01:T;
   x=zeros(6,length(t));
   x(:,1)=x0;
   h=0.01;
   for i=1:length(t)-1
        f1 = RHS_contr(x(:, i), A, B,K);
        f2 = RHS_contr(x(:, i)+f1/2*h, A, B,K);
        f3 = RHS_contr(x(:, i)+f2/2*h, A, B,K);
        f4 = RHS_contr(x(:, i)+f3*h, A, B,K);
        x(:, i+1) = x(:, i) + h*(f1/6+((f2+f3)/3)+f4/6);
   end

    function dx=RHS_contr(x,A,B,K)
        dx = (A-B*K)*x;
   end
end