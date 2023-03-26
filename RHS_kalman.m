function dxdt=RHS_kalman(x,A,L,C,B,K,x_k,x_prim,u_k)
dxdt=(A+L*C)*x + (2*L*C-B*K)*x_k + (L*C+B*K)*x_prim + B*u_k;
end