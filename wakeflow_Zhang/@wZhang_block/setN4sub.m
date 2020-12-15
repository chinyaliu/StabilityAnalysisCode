function [N, arr] = setN4sub(obj)
if isnan(obj.cL)
    obj.cL = obj.zL;
end
if obj.cL > 0.5*obj.zL
    arr = [0 obj.zL 2*obj.zL obj.h];
else
    arr = [0 obj.zL-obj.cL obj.zL obj.zL+obj.cL obj.h];
end
N = zeros(1,length(arr)-1);
% a = diff(arr)/obj.h < 0.2;
% N(a) = ceil(obj.N/5);
a = diff(arr)/obj.h < 0.1;
N(a) = ceil(obj.N/10);
N = N + ceil((obj.N-sum(N))*diff(arr)/obj.h).*(~a);
N(end) = obj.N-sum(N(1:end-1));
end