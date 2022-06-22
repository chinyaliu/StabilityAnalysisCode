function [x_max] = goldensearchmax(xL,xU,f,varargin)
%
% Golden section search method for local maximum of f(x)
% init     : set to TRUE for first iteration
% [xL, xU] : interval to define the local maximum
% f        : function handle f(x,varargin)
% varargin : additional parameters required to define f(x)
%
d = @(xl,xu) 0.5*(sqrt(5)-1)*(xu-xl); % golden ratio
x1 = xL + d(xL,xU);
x2 = xU - d(xL,xU);
err = 1e-5;
f1 = f(x1,varargin{:});
f2 = f(x2,varargin{:});
for i = 1:100
   if (xU-xL) < err
       x_max = 0.5*(xU+xL);
       break;
   end
   if f1 > f2
       xL = x2;
       x2 = x1;
       f2 = f1;
       x1 = xL + d(xL,xU);
       f1 = f(x1,varargin{:});
   else
       xU = x1;
       x1 = x2;
       f1 = f2;
       x2 = xU - d(xL,xU);
       f2 = f(x2,varargin{:});
   end
end
end