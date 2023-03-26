function dXdt=DRE(X,A,B,R,E,Q)
dXdt=-(E^-1)' * A' * X - X * A * E^-1 + X * B * R^-1 * B' * X - (E^-1)' * Q * E^-1;
end