function F = mySystemParasitoid2(X,b,a,lambda,mu)

    x1 = parasitoidMap(X,b,a,lambda,mu);
    x2 = parasitoidMap(x1,b,a,lambda,mu);
  
    F = x2 - X;

end