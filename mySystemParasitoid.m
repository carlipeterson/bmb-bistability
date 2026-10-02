function F = mySystemParasitoid(X,b,a,lambda,mu)

    F = X - parasitoidMap(X,b,a,lambda,mu);

end