function [A, B] = matAB(obj, subD)
ln = length(subD);
A = zeros(obj.N+ln+2,obj.N+ln+1);
B = zeros(obj.N+ln+1,obj.N+ln+1);
switch obj.method(1)
    case 'ray'
       %% Rayleigh
        % Governing equation
        x = 1; y = 1;
        for i = 1:length(subD)
            A(y:y+subD(i).N+2, x:x+subD(i).N) = subD(i).A;
            B(y+2:y+subD(i).N, x:x+subD(i).N) = subD(i).B;
            x = x+subD(i).N+1;
            y = y+subD(i).N+1;
        end
        % BC (free surface)
        [A(1:2,1:subD(1).N+1), A(1:2,end), B(1:2,1:subD(1).N+1), B(1:2,end)] = subD(1).BC0(obj.Fr2);
        % BC (truncated)
        A(end, obj.N-subD(end).N+ln:end-1) = subD(end).BCh();
        A = A(1:end-1,:);
    otherwise
        error('Invalid create AB matrix method name');
end
end