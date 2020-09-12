function [A, B] = matAB_poiseuille(method,N, D, U, BC, k, Re)
% D0 = D(:,:,1); D1 = D(:,:,2); D2 = D(:,:,3); D3 = D(:,:,4); D4 = D(:,:,5);
% U = U(:,:,1); Uz = U(:,:,2); Uzz = U(:,:,3); 
switch lower(method)
    case 'd2'
        %% D2
        % Matrix A
        A11 = 1i*k*U(:,:,3).*D(:,:,1);
        A21 = D(:,:,3)-k^2*D(:,:,1);
        A12 = A21/Re - 1i*k*U(:,:,1).*D(:,:,1);
        A22 = -D(:,:,1);
        BC = [BC{1}(1,:);BC{2}(1,:);BC{1}(2,:);BC{2}(2,:)];
        A = [A11 A12; A21 A22; BC zeros(4,N+1)];
        A = A./(-1i); % eigenvalue= kc
        % Matrix B
        B = zeros(2*(N+1),2*(N+1));
        B(1:N-1,N+2:end) = D(:,:,1);
    case 'uw'
        %% uw
        % Matrix A
        A11 = D(:,:,4)/Re - (k^2/Re + 1i*k*U(:,:,1)).*D(:,:,2) - 1i*k*U(:,:,2).*D(:,:,1);
        A21 = 1i*k*D(:,:,1);
        A12 = -1i*k*D(:,:,3)/Re - U(:,:,2).*D(:,:,2) + (1i*k^3/Re - U(:,:,1)*k^2 - U(:,:,3)).*D(:,:,1);
        A22 = D(:,:,2);
        BC = [BC{1}(1,:);BC{2}(1,:)];
        A = [A11 A12; A21 A22; BC zeros(2,N+1); zeros(2,N+1) BC];
        % Matrix B
        B11 = -1i*D(:,:,2);
        B12 = -k*D(:,:,1);
        B = [B11 B12; zeros(N+3,2*(N+1))];
    case 'd4'
        %% D4
        % Matrix A
        A = 1i*D(:,:,5)/Re + (k*U(:,:,1) - 2*1i*k^2/Re).*D(:,:,3) + (1i*k^4/Re - U(:,:,1)*k^3 - U(:,:,3)*k).*D(:,:,1);
        BC = [BC{1}(1,:);BC{2}(1,:);BC{1}(2,:);BC{2}(2,:)];
        A = [A; BC];
        % Matrix B
        B = D(:,:,3) - k^2*D(:,:,1);
        B = [B; zeros(4,N+1)];
    otherwise
        error('Invalid create AB matrix method name');
end
end