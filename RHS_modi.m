function dxdt=RHS_modi(x,A,B,E,K,u) %solving for x' = the actual trajectori of the dual pendulum vs x_bar = the optimal trajectory
%dxdt=E^-1*(A+B*K)*x;
dxdt=-E^-1 * B * K * x + E^-1 * B * u + E^-1*(A+B*K) * x;
end