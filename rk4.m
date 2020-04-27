function Xnext=rk4(Xnow,dt,params)

k1=deriv(Xnow,params)*dt;
k2=deriv(Xnow+k1/2,params)*dt;
k3=deriv(Xnow+k2/2,params)*dt;
k4=deriv(Xnow+k3,params)*dt;
Xnext=Xnow+(k1+2*k2+2*k3+k4)/6;


