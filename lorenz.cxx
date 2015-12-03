#include <cmath>
#include <iostream>

using namespace std;

void f(double* k, double* k1, double* k2, double* k3, const int a, const int b, const double c, double dt, double* p, double* r){
k[0] = a*((r[1] + dt*p[0]*k[1] + dt*p[1]*k2[1] + dt*p[2]*k3[1]) - (r[0] + dt*p[0]*k1[0] + dt*p[1]*k2[0] + dt*p[2]*k3[0]));
k[1] = (r[0] + dt*p[0]*k1[0] + dt*p[1]*k2[0] + dt*p[2]*k3[0])*(b - (r[2] + dt*p[0]*k1[2] + dt*p[1]*k2[2] + dt*p[2]*k3[2])) - 
	(r[1] + dt*p[0]*k1[1] + dt*p[1]*k2[1] + dt*p[2]*k3[1]);
k[2] = (r[0] + dt*p[0]*k1[0] + dt*p[1]*k2[0] + dt*p[2]*k3[0]) * (r[1] + dt*p[0]*k1[1] + dt*p[1]*k2[1] + dt*p[2]*k3[1]) - c*(r[2] + dt*p[0]*k1[2] + dt*p[1]*k2[2] + dt*p[2]*k3[2]);

}

int main(){
const int N = 100000;
double dt = 100.0/N;
double r[3];
double p[3];
double k1[3];
double k2[3];
double k3[3];
double k4[3];
int a = 10;
int b = 28;
double c = 8/3;
r[0]=1.0;
r[1]=1.0;
r[2]=1.0;

//r1[0] = a*(y-x);
//r1[1] = x*(b-z) - y;
//r1[2] = x*y - c*z;

for(int i=0; i<N; i++){
p[0]=0; p[1]=0; p[2]=0;
f(k1, k1, k2, k3, a, b, c, dt, p, r);
p[0]=0.5; p[1]=0; p[2]=0;
f(k2, k1, k2, k3, a, b, c, dt, p, r);
p[0]=0; p[1]=0.5; p[2]=0;
f(k3, k1, k2, k3, a, b, c, dt, p, r);
p[0]=0; p[1]=0; p[2]=1;
f(k4, k1, k2, k3, a, b, c, dt, p, r);

r[0] = r[0] + dt/6*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0]);
r[1] = r[1] + dt/6*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1]);
r[2] = r[2] + dt/6*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2]);

cout << i*dt << "\t" << r[0] << "\t" << r[1] << "\t" << r[2] << endl;
}
return 0;
}
