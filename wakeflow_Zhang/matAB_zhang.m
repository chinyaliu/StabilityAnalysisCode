function [A, B] = matAB_zhang(method,N, D, U, BC, k, Re, Fr2)
% D0 = D(:,:,1); D1 = D(:,:,2); D2 = D(:,:,3); D3 = D(:,:,4); D4 = D(:,:,5);
% U = U(:,:,1); Uz = U(:,:,2); Uzz = U(:,:,3); 
switch lower(method)
    case 'd2'
    %% D2
        % Matrix A
        A11 = (1i*k^4/Re - U(:,:,1)*k^3 - U(:,:,3)*k).*D(:,:,1);
        A12 = 1i*D(:,:,3)/Re + (k*U(:,:,1) - 2*1i*k^2/Re).*D(:,:,1);
        A21 = D(:,:,3); A22 = -D(:,:,1);
        Abc1 = [k*BC{1}(1,:) zeros(1,N+1) k*BC{3}(1,1)];
        Abc2 = [k^2*BC{1}(1,:) BC{1}(1,:) BC{3}(3,1)];
        Abc3 = [(k*BC{3}(1,:)-3*1i*k^2/Re).*BC{1}(2,:)-k*BC{3}(2,:).*BC{1}(1,:) 1i*BC{1}(2,:)/Re k/Fr2];
        A = [[A11 A12; A21 A22] zeros(2*N-2,1);...
            Abc1; Abc2; Abc3; [BC{2}(1,:) zeros(1,N+1) 0]; [zeros(1,N+1) BC{2}(1,:) 0]];
        % Matrix B
        B = zeros(2*N+3, 2*N+3);
        B(1:N-1,1:N+1) = -k^2*D(:,:,1);
        B(1:N-1,N+2:end-1) = D(:,:,1);
        B(end-2,1:N+1) = BC{1}(2,:);
        B(2*N-1,end) = 1;
    case 'uw'
        %% uw
        % Matrix A
        A11 = -D(:,:,4)/Re + (k^2/Re + 1i*k*U(:,:,1)).*D(:,:,2) + 1i*k*U(:,:,2).*D(:,:,1);
        A12 = 1i*k*D(:,:,3)/Re + U(:,:,2).*D(:,:,2) + (U(:,:,3) + U(:,:,1)*k^2 - 1i*k^3/Re).*D(:,:,1);
        A21 = 1i*k*D(:,:,1); A22 = D(:,:,2);
        Abc1 = [zeros(1,N+1) 1i*BC{1}(1,:) k*BC{3}(1,1)];
        Abc2 = [BC{1}(2,:) 1i*k*BC{1}(1,:) BC{3}(3,1)];
        Abc3 = [1i*BC{1}(3,:)/Re+(k*BC{3}(1,:)-1i*k^2/Re).*BC{1}(1,:) 2*k*BC{1}(2,:)/Re-BC{3}(2,:).*BC{1}(1,:) k/Fr2];
        A = [[A11 A12; A21 A22] zeros(2*N-2,1);Abc1; Abc2; Abc3; ...
            [BC{2}(2,:) zeros(1,N+2)]; [zeros(1,N+1) BC{2}(1,:) 0]];
        % Matrix B
        B = zeros(2*N+3, 2*N+3);
        B(1:N-1,1:N+1) = 1i*D(:,:,2);
        B(1:N-1,N+2:end-1) = k*D(:,:,1);
        B(end-2,1:N+1) = BC{1}(1,:);
        B(2*N-1,end) = 1;
    case 'd4'
      %% D4
        % Matrix A
        A1 = 1i*D(:,:,5)/Re + (k*U(:,:,1) - 2*1i*k^2/Re).*D(:,:,3) + (1i*k^4/Re - U(:,:,1)*k^3 - U(:,:,3)*k).*D(:,:,1);
        Abc1 = [k*BC{1}(1,:) k*BC{3}(1,1)];
        Abc2 = [BC{1}(3,:) + k^2*BC{1}(1,:) BC{3}(3,1)];
        Abc3 = [1i*BC{1}(4,:)/Re + (k*BC{3}(1,:) - 3*1i*k^2/Re).*BC{1}(2,:) - k*BC{3}(2,:).*BC{1}(1,:) k/Fr2];
        A = [[A1 zeros(N-3,1)]; Abc1; Abc2; Abc3; [BC{2}(1,:) 0]; [BC{2}(3,:) 0]];
        % Matrix B
        B = zeros(N+2, N+2);
        B(1:N-3,1:N+1) = D(:,:,3) - k^2*D(:,:,1);
        B(N,1:N+1) = BC{1}(2,:);
        B(N-2,end) = 1;
    otherwise
        error('Invalid create AB matrix method name');
end
end