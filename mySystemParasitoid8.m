function F = mySystemParasitoid8(X,b,a,lambda,mu)

    x1 = parasitoidMap(X,b,a,lambda,mu);
    x2 = parasitoidMap(x1,b,a,lambda,mu);
    x3 = parasitoidMap(x2,b,a,lambda,mu);
    x4 = parasitoidMap(x3,b,a,lambda,mu);
    x5 = parasitoidMap(x4,b,a,lambda,mu);
    x6 = parasitoidMap(x5,b,a,lambda,mu);
    x7 = parasitoidMap(x6,b,a,lambda,mu);
    x8 = parasitoidMap(x7,b,a,lambda,mu);
  
    F = x8 - X;

end