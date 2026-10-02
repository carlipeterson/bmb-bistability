function J_fun = makeJacobianParasitoid(b,a,lambda,mu)
% makeJacobian  Create a numeric function handle for the Jacobian
%   J_fun = makeJacobian(params)

    pop_syms = sym('x',[2 1],'real');
    
    
    x = pop_syms(1);
    y = pop_syms(2);
    
    % Define the system

    % F1 = (r+1)*x-r*x^2-c*x*y;
    % F2 = c*x*y;

    % F1 = (r+1)*x-r*x^2-c*x*y;
    % F2 = c*x*y+(1-d)*y;

    % F1 = x*exp(r*(1-x-y));
    % F2 = c*x*y;

    % F1 = x*exp(r*(1-x)-y);
    % F2 = c*x*(1-exp(-y));

    % F1 = (r+1)*x-r*x^2-c*x*y/(x+1);
    % F2 = c*x*y/(x+1);

    % F1 = x*exp(r-x-y);
    % F2 = c*x*(1-exp(-y));

    % F1 = lambda*x*exp(-mu*x)*exp(-a*y);
    % F2 = b*x*(1-exp(-a*y));

    F1 = lambda*x*exp(-mu*x)*exp(-a*y);
    F2 = b*lambda*x*exp(-mu*x)*(1-exp(-a*y));

    
    F = [F1; F2];

    all_vars = [pop_syms];
    
    J = jacobian(F,all_vars);
    
    J_fun = matlabFunction(J,'Vars',{all_vars});

end
