function [A, B] = matAB(obj, U)
% D0 = D(:,:,1); D1 = D(:,:,2); D2 = D(:,:,3); D3 = D(:,:,4); D4 = D(:,:,5);
% U = U(:,1); Uz = U(:,2); Uzz = U(:,3); 
N = obj.N; k = obj.k; D = obj.D;
A = zeros(N+2,N+2);
B = zeros(N+2,N+2);
switch obj.method(1)
    case 'ray'
    %% Rayleigh
        % Governing equation
        A_ge = k*U(:,1).*D(:,:,3) - (U(:,1)*k^3 + U(:,3)*k).*D(:,:,1);
        A(1:N-1, 1:N+1) = A_ge(2:N,:);
        B_ge = D(:,:,3) - k^2*D(:,:,1);
        B(1:N-1, 1:N+1) = B_ge(2:N,:);
        % KFSBC
        A(end-2, 1:N+1) = k*D(1,:,1);
        A(end-2, end) = k*U(1,1);
        B(end-2, end) = 1;
        % DFSBC
        A(end-1, 1:N+1) = k*U(1,1).*D(1,:,2) - k*U(1,2).*D(1,:,1);
        A(end-1, end) = k/obj.Fr2;
        B(end-1, 1:N+1) = D(1,:,2);
        % z = -h
        A(end, 1:N+1) = D(end,:,2)-k*D(end,:,1);
%     case 'd4'
%       %% D4
%         % Governing equation
%         A_ge = 1i*D(:,:,5)/obj.Re + (k*U(:,1) - 2*1i*k^2/obj.Re).*D(:,:,3) + (1i*k^4/obj.Re - U(:,1)*k^3 - U(:,3)*k).*D(:,:,1);
%         B_ge = D(:,:,3) - k^2*D(:,:,1);
%         switch obj.method(3)
%             case 'schimd'
%                 A(1:N-3, 1:N+1) = A_ge(3:N-1,:);
%                 A(N-2:2*N-6, N+2:end-1) = A_ge(N+4:end-2,:);
%                 B(1:N-3, 1:N+1) = B_ge(3:N-1,:);
%                 B(N-2:2*N-6, N+2:end-1) = B_ge(N+4:end-2,:);
%                 % Matching condition
%                 A(end-8:end-5, 1:end-1) = permute([D(N+1,:,1:4) -D(N+2,:,1:4)],[3 2 1]);
%             case 'd4'
%                 A(1:N-3, 1:N+1) = A_ge(2:N-2,:);
%                 A(N-2:2*N-6, N+2:end-1) = A_ge(N+1:end-1,:);
%                 B(1:N-3, 1:N+1) = B_ge(2:N-2,:);
%                 B(N-2:2*N-6, N+2:end-1) = B_ge(N+1:end-1,:);
%                 % Matching condition
%                 A(end-8:end-5, 1:end-1) = permute([D(N-1,:,1:4) -D(N,:,1:4)],[3 2 1]);
%         end
%         % KFSBC
%         A(end-4, 1:N+1) = k*D(1,:,1);
%         A(end-4, end) = k*U(1,1);
%         B(end-4, end) = 1;
%         % DFSBC tangent
%         A(end-3, 1:N+1) = D(1,:,3) + k^2*D(1,:,1);
%         A(end-3, end) = U(1,3);
%         % DFSBC normal
%         A(end-2, 1:N+1) = 1i*D(1,:,4)/obj.Re + (k*U(1,1) - 3*1i*k^2/obj.Re).*D(1,:,2) - k*U(1,2).*D(1,:,1);
%         A(end-2, end) = k/obj.Fr2;
%         B(end-2, 1:N+1) = D(1,:,2);
%         % z = -h
%         A(end-1, N+2:end-1) = D(end,:,1);
%         A(end, N+2:end-1) = D(end,:,3);
    otherwise
        error('Invalid create AB matrix method name');
end
end