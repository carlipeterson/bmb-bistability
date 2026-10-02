function xnext = parasitoidMap(X,b,a,lambda,mu)

    x = X(1);
    y = X(2);

    xnext = zeros(2,1);

    % xnext(1) = (r+1)*x-r*x^2-c*x*y;
    % xnext(2) = c*x*y;

    % xnext(1) = (r+1)*x-r*x^2-c*x*y;
    % xnext(2) = c*x*y+(1-d)*y;

    % xnext(1) = x*exp(r*(1-x-y));
    % xnext(2) = c*x*y;

    % xnext(1) = x*exp(r*(1-x)-y);
    % xnext(2) = c*x*(1-exp(-y));

    % xnext(1) = (r+1)*x-r*x^2-c*x*y/(x+1);
    % xnext(2) = c*x*y/(x+1);

    % xnext(1) = x*exp(r-x-y);
    % xnext(2) = c*x*(1-exp(-y));

    % xnext(1) = lambda*x*exp(-mu*x)*exp(-a*y);
    % xnext(2) = b*x*(1-exp(-a*y));

    xnext(1) = lambda*x*exp(-mu*x)*exp(-a*y);
    xnext(2) = b*lambda*x*exp(-mu*x)*(1-exp(-a*y));
    

end