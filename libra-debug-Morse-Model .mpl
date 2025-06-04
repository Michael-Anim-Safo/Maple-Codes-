
restart; Digits:=16;
ix:=x=-infinity..infinity;

#V:=160*(exp(-1.02*(x-1.4))-1)^2;  #Morse
V := 160*(1 - exp(-1.02*(x - 1.4)))^2



;
plot(V,x=0.5..5);
a:=9.16*7;
q1:= 0.9454764334708259;
q2:= 1.6545235665291742 ;
g1:=exp(-a*(x-q1)^2/2)*(a/Pi)^(1/4);
g2:=exp(-a*(x-q2)^2/2)*(a/Pi)^(1/4);
K11:=int(diff(g1,x)^2/2,x=-infinity..infinity);K12:=int(diff(g1,x)*diff(g2,x)/2,ix);
K21:=int(diff(g1,x)*diff(g2,x)/2,x=-infinity..infinity);K22:=int(diff(g2,x)*diff(g2,x)/2,ix);
qc:=(q1+q2)/2;
#Vh:=convert(series(V,x=qc,3),polynom):
V1:=convert(series(V,x=q1,3),polynom);
V2:=convert(series(V,x=q2,3),polynom);
Vh:=(V1+V2)/2:


print('EXACT_MORSE_INTERGRALS');
v11:=int(g1^2*V,ix)+K11;v12:=int(g1*g2*V,ix)+K12;
v21:=int(g1*g2*V,ix)+K21;v22:=int(g2*g2*V,ix)+K22;


print('AVERAGE_of_LHAs');
v11:=int(g1^2*V1,ix)+K11;  v12:=int(g1*g2*Vh,ix)+K12;
v21:=int(g1*g2*Vh,ix)+K21; v22:=int(g2*g2*V2,ix)+K22;

Vh:=convert(series(V,x=qc,2),polynom):
V1:=convert(series(V,x=q1,2),polynom);
V2:=convert(series(V,x=q2,2),polynom);
Vh:=(V1+V2)/2:

print('BAT');
v11:=int(g1^2*V1,ix)+K11;  v12:=int(g1*g2*Vh,ix)+K12;
v21:=int(g1*g2*Vh,ix)+K21; v22:=int(g2*g2*V2,ix)+K22;


plot({V1,V2,V},x=0.5..5,y=-10..200);



