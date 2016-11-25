function [x,i,v] = newton_equ(A, b, x_1 , f, df, Hession)
%% 
%  find the smallest value of f(x),
%  with constraint Ax = b, and initialization x_1.
%  by Newton method

%  df  is the gradient of f(x),  and Hesstion is the Hesstion of f(x)
% 

 RXV = @(x,v) [df(x) + A'*v ; A*x - b];

%%

[m, n] = size(A);

maxIter = 200;

x = x_1;

theta = 1.0e-10;

v = zeros(m,1);
alpha = 0.4;
beta = 0.7;
for i = 1:maxIter
    
    r = [x;v]; 
    grad = [df(x) + A'*v; A*x - b];
    H = [Hession(x),A';A,zeros(m)];
    d_newton = -H\grad;
    if (norm(grad) < theta )
        disp('over')
        
        break;
    end
    tempx = d_newton(1:n);
    max_step = x./(abs(tempx) + theta/100);

    t = min(max_step(tempx < 0))/1.1;
    t = min(1,t);
    R0 = RXV( x , v ) ;
    R0 = (R0'*R0)^0.5;
    tempR = RXV( x + t*d_newton(1:n), v + t*d_newton(n+1:end));
    
    tempR = (tempR'*tempR)^0.5;
    while(tempR> (1 - alpha*t)*R0)
        t = t*beta;
     tempR = RXV( x + t*d_newton(1:n), v + t*d_newton(n+1:end));
        
        tempR = (tempR'*tempR)^0.5;
    end
    
    r = r+t*d_newton;
    x = r(1:n);
    v = r(n+1:end);    

end
