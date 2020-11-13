function [A, B] = matAB_zhang2(obj, D, U, w1, w2)
% D0 = D(:,:,1); D1 = D(:,:,2); D2 = D(:,:,3); D3 = D(:,:,4); D4 = D(:,:,5);
% U = U(:,1); Uz = U(:,2); Uzz = U(:,3); 
N = obj.N; k = obj.k;
switch lower(obj.method(1))
    case 'ray'
    %% Rayleigh
        % Matrix A
        A1 = [k*U{1}(:,1).*D{1}(:,:,3) - (U{1}(:,1)*k^3 + U{1}(:,3)*k).*D{1}(:,:,1) zeros(N-1,N+2)];
        A2 = [zeros(N-1,N+1) k*U{2}(:,1).*D{2}(:,:,3) - (U{2}(:,1)*k^3 + U{2}(:,3)*k).*D{2}(:,:,1) zeros(N-1,1)];
        Abc1 = [k*obj.BC{1}(1,:) zeros(1,N+1) k*obj.BC{3}(1)];
        Abc3 = [k*obj.BC{3}(1).*obj.BC{1}(2,:)*w1(2) - k*obj.BC{3}(2).*obj.BC{1}(1,:) zeros(1,N+1) k/obj.Fr2];
        Amatch = [obj.BC{2}(1:2,:).*w1(1:2).' -obj.BC{1}(1:2,:).*w2(1:2).' zeros(2,1)];
        A = [A1; A2; Abc1; Abc3; [zeros(1,N+1) obj.BC{2}(1,:) 0]; Amatch];
        % Matrix B
        B = zeros(2*N+3, 2*N+3);
        B(1:N-1,1:N+1) = D{1}(:,:,3) - obj.k^2*D{1}(:,:,1);
        B(N:2*N-2,N+2:end-1) = D{2}(:,:,3) - k^2*D{2}(:,:,1);
        B(2*N,1:N+1) = w1(2)*obj.BC{1}(2,:);
        B(2*N-1,end) = 1;
    case 'd4'
      %% D4
        % Matrix A
        A1 = [1i*D{1}(:,:,5)/obj.Re + (k*U{1}(:,1) - 2*1i*k^2/obj.Re).*D{1}(:,:,3) + (1i*k^4/obj.Re - U{1}(:,1)*k^3 - U{1}(:,3)*k).*D{1}(:,:,1) zeros(N-3,N+2)];
        A2 = [zeros(N-3,N+1) 1i*D{2}(:,:,5)/obj.Re + (k*U{2}(:,1) - 2*1i*k^2/obj.Re).*D{2}(:,:,3) + (1i*k^4/obj.Re - U{2}(:,1)*k^3 - U{2}(:,3)*k).*D{2}(:,:,1) zeros(N-3,1)];
        Abc1 = [k*w1(1)*obj.BC{1}(1,:) zeros(1,N+1) k*obj.BC{3}(1,1)];
        Abc2 = [w1(3)*obj.BC{1}(3,:) + k^2*w1(1)*obj.BC{1}(1,:) zeros(1,N+1) obj.BC{3}(3)];
        Abc3 = [1i*w1(4)*obj.BC{1}(4,:)/obj.Re + (k*obj.BC{3}(1) - 3*1i*k^2/obj.Re).*obj.BC{1}(2,:)*w1(2) - k*obj.BC{3}(2).*obj.BC{1}(1,:)*w1(1) zeros(1,N+1) k/obj.Fr2];
        Amatch = [obj.BC{2}(1:4,:).*w1(1:4).' -obj.BC{1}(1:4,:).*w2(1:4).' zeros(4,1)];
        A = [A1; A2; Abc1; Abc2; Abc3; [zeros(1,N+1) w2(1)*obj.BC{2}(1,:) 0]; [zeros(1,N+1) w2(3)*obj.BC{2}(3,:) 0]; Amatch];
        % Matrix B
        B = zeros(2*N+3, 2*N+3);
        B(1:N-3,1:N+1) = D{1}(:,:,3) - k^2*D{1}(:,:,1);
        B(N-2:2*N-6,N+2:end-1) = D{2}(:,:,3) - k^2*D{2}(:,:,1);
        B(2*N-3,1:N+1) = w1(2)*obj.BC{1}(2,:);
        B(2*N-5,end) = 1;
    otherwise
        error('Invalid create AB matrix method name');
end
end