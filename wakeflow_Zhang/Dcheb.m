function [x,D] = Dcheb(N,num,points)
if strcmpi(points,'d4') M = N-2; else M = N; end
x = cospi((0:1:M)/M)'; 
D = zeros(M+1,N+1,num+1);
D(:,:,1) = cospi((0:1:M)'*(0:1:N)/M);
D(:,2:3,2) = [D(:,1,1) 4*D(:,2,1)];
D(:,3,3) = 4*D(:,1,1);
for i = 4:N+1
    D(:,i,2:end) = (i-1)/(i-3)*D(:,i-2,2:end) + 2*(i-1)*D(:,i-1,1:end-1);
end
end