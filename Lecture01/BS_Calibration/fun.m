function dp=fun(vol,spot,strike,rf,maturity,pmkt)
%compute B&S Price
d1=(log(spot./strike)+(rf+vol.*vol/2).*maturity)./...
    (vol*sqrt(maturity));
d2=d1-vol.*sqrt(maturity);
nd1=normcdf(d1,0,1); nd2=normcdf(d2,0,1);
bs=(spot.*nd1-strike.*exp(-rf.*maturity).*nd2);
%compute the difference between the
%B&S price and the real price
dp=bs-pmkt;