$THETA
0.673 14.3 3.05

$CMT GUT CENT


$MAIN
double KA = THETA1*exp(ETA(1));
double VC = THETA2*exp(ETA(2));
double CL = THETA3*exp(ETA(3));

$ODE
dxdt_GUT = -KA*GUT;
dxdt_CENT = KA*GUT - (CL/VC)*CENT;


$TABLE
capture IPRED = CENT/VC;
double DV = IPRED * (1 + EPS(1)) + EPS(2);

$CAPTURE DV


$OMEGA
0.1
0.0958
0.0961

$SIGMA
0.106 
0

$SET delta=0.1, end=24*3