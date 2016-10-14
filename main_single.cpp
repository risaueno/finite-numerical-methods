#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

double energy(double,double);
void print(double,double,double,ofstream&);
void stability(double,double,double,ofstream&);
void euler(double,double,double,double,int,ofstream&);
void leapfrog(double,double,double,double,int,ofstream&);
void rungekutta(double,double,double,double,int,ofstream&);
void impliciteuler(double,double,double,double,int,ofstream&);
void eulerS(double,double,double,double,int,ofstream&);
void leapfrogS(double,double,double,double,int,ofstream&);
void rungekuttaS(double,double,double,double,int,ofstream&);
void impliciteulerS(double,double,double,double,int,ofstream&);
void criticaldt(double,double,double,ofstream&,ofstream&,ofstream&,ofstream&);

/*	MAIN
--------------------------------------------------------------------------------*/
int main()
{
    ofstream stream01("Single Undamped Euler.txt");
    ofstream stream02("Single Undamped Leapfrog.txt");
    ofstream stream03("Single Undamped Runge-Kutta.txt");
    ofstream stream04("Single Undamped Implicit Euler.txt");
    ofstream stream05("Single Damped Euler.txt");
    ofstream stream06("Single Damped Leapfrog.txt");
    ofstream stream07("Single Damped Runge-Kutta.txt");
    ofstream stream08("Single Damped Implicit Euler.txt");

    ofstream stream01S("Stability Single Undamped Euler.txt");
    ofstream stream02S("Stability Single Undamped Leapfrog.txt");
    ofstream stream03S("Stability Single Undamped Runge-Kutta.txt");
    ofstream stream04S("Stability Single Undamped Implicit Euler.txt");
    ofstream stream05S("Stability Single Damped Euler.txt");
    ofstream stream06S("Stability Single Damped Leapfrog.txt");
    ofstream stream07S("Stability Single Damped Runge-Kutta.txt");
    ofstream stream08S("Stability Single Damped Implicit Euler.txt");

	//Undamped
    double omega=0, theta=0.1, D=0, dt=0.1;
    int N=1000/dt;
    euler(omega,theta,D,dt,N,stream01);   //Prints info for each method
    leapfrog(omega,theta,D,dt,N,stream02);
    rungekutta(omega,theta,D,dt,N,stream03);
    impliciteuler(omega,theta,D,dt,N,stream04);
	criticaldt(D,omega,theta,stream01S,stream02S,stream03S,stream04S);  //Outputs stability check file for all methods

	//Damped
    D=0.2;
    euler(omega,theta,D,dt,N,stream05);
    leapfrog(omega,theta,D,dt,N,stream06);
    rungekutta(omega,theta,D,dt,N,stream07);
    impliciteuler(omega,theta,D,dt,N,stream08);
	criticaldt(D,omega,theta,stream05S,stream06S,stream07S,stream08S);

    return 0;
}

/*	FUNCTIONS
--------------------------------------------------------------------------------*/

//Energy calculation
double energy(double theta, double omega)
{
    return theta*theta + omega*omega;
}

//Print t, theta and total energy
void print(double t, double theta, double omega, ofstream& stream)
{
    stream << t << ' ' << theta << ' ' << energy(theta, omega) << endl;
}

//Outputs table for stability of dt
void stability(double Ei, double Ef, double dt, ofstream& stream)
{
    if(Ef > 10000000*Ei)  //Checks if maximum energy at a later time is  significantly greater than the initial energy
    { stream << dt << ' ' << 1 << endl; } //Output 1 if unstable
    else
    { stream << dt << ' '<< 0 << endl; } //Output 0 if stable
}

//Explicit Euler
void euler(double omega, double theta,  double D, double dt, int N, ofstream& stream)
{
    double t=0, omega0=0;
    for(int i=0; i<N ; i++)
    {
        print(t, theta, omega, stream); //Output of time, theta and energy
        omega0= omega;
        omega-= (theta + D*omega)*dt;
        theta+= omega0*dt;
        t+=dt;
    }
}

//Leapfrog
void leapfrog(double omega, double theta,  double D, double dt, int N, ofstream& stream)
{
    double t=0, omega0=omega, theta0=theta, omega1=0, theta1=0;
    print(t, theta, omega, stream);
    omega-= (theta + D*omega)*dt;   //First step: Euler
    theta+= omega0*dt;

    for (int i=1; i<N; i++)    //Subsequent steps (loop): Leapfrog
    {
        omega1=omega, theta1=theta;
        t+=dt;
        print(t, theta, omega, stream);
        omega=omega0-2*(theta+D*omega)*dt; //Update omega and theta
        theta=theta0+2*omega1*dt;
        theta0=theta1, omega0=omega1; //Update initial values
    }
}

//Runge-Kutta
void rungekutta(double omega, double theta,  double D, double dt, int N, ofstream& stream)
{
    double t=0, p1=0, p2=0, p3=0, p4=0, k1=0, k2=0, k3=0, k4=0;
    for(int i=1; i<N; i++)
    {
        print(t, theta, omega, stream);
        p1= dt*omega;
        k1= -dt*(theta+D*omega);
        p2= dt*(omega+0.5*k1);
        k2= -dt*(theta+0.5*p1+D*(omega+0.5*k1));
        p3= dt*(omega+0.5*k2);
        k3= -dt*(theta+0.5*p2+D*(omega+0.5*k2));
        p4= dt*(omega+k3);
        k4= -dt*(theta+p3+D*(omega+k3));
        theta= theta+(p1+2*p2+2*p3+p4)/6;
        omega= omega+(k1+2*k2+2*k3+k4)/6;
        t+=dt;
    }
}

//Implicit Euler
void impliciteuler(double omega, double theta,  double D, double dt, int N, ofstream& stream)
{
    double t=0, a=0, b=0;
    for(int i=1; i<N; i++)
    {
        print(t, theta, omega, stream);
        a= 1+D*dt, b= a+dt*dt;
        theta= theta*a/b + omega*dt/b;  //Calculates new theta
        omega= omega/a - dt*theta/a;    //Omega depends on new theta
        t+=dt;
    }
}

//Euler stability test
void eulerS(double omega, double theta,  double D, double dt, int N, ofstream& stream)
{
    double omega0=0;
    double Ei= energy(theta, omega), Ef=0;  //Ei = initial energy, Ef = final energy
    for(int i=0; i<N ; i++)
    {
        omega0= omega;
        omega-= (theta + D*omega)*dt;
        theta+= omega0*dt;
        if (energy(theta, omega)>Ef)    //Loop maintains maximum value of Ef
        {
            Ef = energy(theta, omega);
        }
    }
    stability(Ei,Ef,dt,stream);
}

//Leapfrog stability test
void leapfrogS(double omega, double theta,  double D, double dt, int N, ofstream& stream)
{
    double omega0=omega, theta0=theta, omega1=0, theta1=0;
    double Ei= energy(theta, omega), Ef=0;
    omega-= (theta + D*omega)*dt;
    theta+= omega0*dt;
    for (int i=1; i<N; i++)
    {
        omega1=omega, theta1=theta;
        omega=omega0-2*(theta+D*omega)*dt;
        theta=theta0+2*omega1*dt;
        theta0=theta1, omega0=omega1;
        if (energy(theta, omega)>Ef)
        {
            Ef = energy(theta, omega);
        }
    }
    stability(Ei,Ef,dt,stream);
}

//Runge-Kutta stability test
void rungekuttaS(double omega, double theta,  double D, double dt, int N, ofstream& stream)
{
    double p1=0, p2=0, p3=0, p4=0, k1=0, k2=0, k3=0, k4=0;
    double Ei= energy(theta, omega), Ef=0;
    for(int i=1; i<N; i++)
    {
        p1= dt*omega;
        k1= -dt*(theta+D*omega);
        p2= dt*(omega+0.5*k1);
        k2= -dt*(theta+0.5*p1+D*(omega+0.5*k1));
        p3= dt*(omega+0.5*k2);
        k3= -dt*(theta+0.5*p2+D*(omega+0.5*k2));
        p4= dt*(omega+k3);
        k4= -dt*(theta+p3+D*(omega+k3));
        theta= theta+(p1+2*p2+2*p3+p4)/6;
        omega= omega+(k1+2*k2+2*k3+k4)/6;
        if (energy(theta, omega)>Ef)
        {
            Ef = energy(theta, omega);
        }
    }
    stability(Ei,Ef,dt,stream);
}

//Implicit Euler stability test
void impliciteulerS(double omega, double theta,  double D, double dt, int N, ofstream& stream)
{
    double a=0, b=0;
    double Ei= energy(theta, omega), Ef=0;
    for(int i=1; i<N; i++)
    {
        a= 1+D*dt, b= a+dt*dt;
        theta= theta*a/b + omega*dt/b;
        omega= omega/a - dt*theta/a;
        if (energy(theta, omega)>Ef)
        {
            Ef = energy(theta, omega);
        }
    }
    stability(Ei,Ef,dt,stream);
}

//Stability function loop
void criticaldt(double D, double omega, double theta, ofstream& stream1, ofstream& stream2, ofstream& stream3, ofstream& stream4)
{
	int N=0;
    for(double dt=0.001; dt <=5; dt+=0.001) //Tests each dt for all methods
    {
        N= ceil(1000/dt); //Runs for time=1000
        eulerS(omega,theta,D,dt,N,stream1);
        leapfrogS(omega,theta,D,dt,N,stream2);
        rungekuttaS(omega,theta,D,dt,N,stream3);
        impliciteulerS(omega,theta,D,dt,N,stream4);
    }
}
