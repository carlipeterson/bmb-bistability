function F = mySystemParasitoid4(X,b,a,lambda,mu)

    x1 = parasitoidMap(X,b,a,lambda,mu);
    x2 = parasitoidMap(x1,b,a,lambda,mu);
    x3 = parasitoidMap(x2,b,a,lambda,mu);
    x4 = parasitoidMap(x3,b,a,lambda,mu);
  
    F = x4 - X;

end