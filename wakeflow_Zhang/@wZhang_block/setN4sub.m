function [N, arr] = setN4sub(obj)
if obj.eps > 0.5*obj.zL
    arr = [0 obj.zL obj.zL+obj.eps obj.h];
    obj.eps = obj.zL;
else
    arr = [0 obj.zL-obj.eps obj.zL obj.zL+obj.eps obj.h];
end
N = zeros(1,length(arr)-1);
a = diff(arr)/obj.h < 0.1;
N(a) = ceil(obj.N/10);
N = N + ceil((obj.N-sum(N))*diff(arr)/obj.h).*(~a);
N(end) = obj.N-sum(N(1:end-1));
% if eps/obj.h < 0.05
%     N(2) = ceil(obj.N/20);
% else
%     N(2) = ceil(obj.eps/obj.h);
% end
% N(1) = ceil((obj.N - 2*N(2))*(obj.zL - obj.eps)/(obj.h-2*obj.eps));
% N(3) = N(2);
% N(4) = obj.N - N(1) - 2*N(2);
% if N(1) < ceil(obj.N/20)
%     N(1) = 0;
%     N(4) = obj.N - 2*N(2);
%     obj.eps = obj.zL;
% end
end