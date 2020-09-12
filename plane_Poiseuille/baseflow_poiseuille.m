function U = baseflow_poiseuille(z,N)
U(:,:,1) = (1-z.^2)*ones(1,N+1);
U(:,:,2) = -2*z*ones(1,N+1);
U(:,:,3) = -2*ones(length(z),N+1);
end