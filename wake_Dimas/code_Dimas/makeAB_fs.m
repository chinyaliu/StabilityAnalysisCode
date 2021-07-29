function [A,B] = makeAB_fs(k,N,h,Fr)
    z = linspace(0,-h,N).';
    [U,Uzz] = baseflow(z);
    dz = z(2)-z(1);
    % A
    a1 = -k*(2*U + Uzz*dz^2 + k*k*U*dz^2);
    a1(end) = a1(end) + U(end)*k*exp(k*dz);
    a2 = U*k;
    a3 = a2;
    agov = spdiags([a2 a1 a3],0:1:2,N,N+3);
    abc_k = [0 k zeros(1,N) a2(1)];
    ad = 0.5*Fr^2*k*U(1)/dz;
    abc_d = [-ad 0 ad zeros(1,N-1) k];
    abc_h = [zeros(1,N) 1 0 0];
    A = full([agov; abc_k; abc_d; abc_h]);
    % B
    no = ones(N,1);
    b1 = (-2-(k*dz)^2)*no;
    b1(end) = b1(end) + exp(k*dz);
    bgov = spdiags([no,b1,[ones(N-1,1);0]],0:1:2,N,N+3);
    bbc_k = [zeros(1,N+2),1];
    bd = 0.5*Fr^2/dz;
    bbc_d = [-bd 0 bd zeros(1,N)];
    bbc_h = zeros(1,N+3);
    B = full([bgov; bbc_k; bbc_d; bbc_h]);
    
    function [U,Uzz] = baseflow(z)
        c1 = 0.9988; c2 = 0.8814;
        U = (1-c1*cosh(c2*z).^(-2));
        Uzz = 2*c1*c2^2*(sech(c2*z)).^2.*((sech(c2*z)).^2-2*(tanh(c2*z)).^2);
    end
end