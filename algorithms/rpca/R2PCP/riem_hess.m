function Hd=riem_hess(d,g,U,V,S)

v1=g-U*(U'*g);
Hd=d+(v1*((d'-V*(V'*d'))*U))*(S\V')...
    +U*(S\((V'*d')*(v1-(v1*V)*V')));

return
