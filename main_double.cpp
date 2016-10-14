#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

double energy(double,double,double,double,double);
void stability(double,double,double,ofstream&);
void print(double,double,double,double,double,double,ofstream&);
double ftheta(double);
double fphi(double);
double fomega(double,double,double,double,double);
double fnu(double,double,double,double,double,double);
void rungekutta(double,double,int,double,double,double,double,double,ofstream&);
void rungekuttaS(double,double,int,double,double,double,double,double,double,ofstream&);
void criticaldt(double,double,double,double,double,double,ofstream&);


/*	MAIN
--------------------------------------------------------------------------------*/

int main()
{
    ofstream stream1("Double Undamped R=1.txt");
    ofstream stream2("Double Undamped R=0.01.txt");
    ofstream stream3("Double Undamped R=100.txt");
    ofstream stream4("Double Damped R=1.txt");
    ofstream stream5("Double Damped R=0.01.txt");
    ofstream stream6("Double Damped R=100.txt");
    ofstream streamSU1("Double Undamped Stability R=1.txt");
    ofstream streamSU2("Double Undamped Stability R=0.01.txt");
    ofstream streamSU3("Double Undamped Stability R=100.txt");
    ofstream streamSD1("Double Damped Stability R=1.txt");
    ofstream streamSD2("Double Damped Stability R=0.01.txt");
    ofstream streamSD3("Double Damped Stability R=100.txt");

    //Undamped
    int N=150000;
    double omega=0, theta=0.1, nu=0, phi=0, G=0, dt=0.001;

    rungekutta(1,G,N,dt,theta,phi,omega,nu,stream1);            //Output info
    rungekutta(0.01,G,N,dt,theta,phi,omega,nu,stream2);
    rungekutta(100,G,N,dt,theta,phi,omega,nu,stream3);
    criticaldt(1,G,theta,phi,omega,nu,streamSU1);               //Stability test for critical timestep
    criticaldt(0.01,G,theta,phi,omega,nu,streamSU2);
    criticaldt(100,G,theta,phi,omega,nu,streamSU3);

    //Damped
    G=1;
    rungekutta(1,G,N,dt,theta,phi,omega,nu,stream4);            //Output info
    rungekutta(0.01,G,N,dt,theta,phi,omega,nu,stream5);
    rungekutta(100,G,N,dt,theta,phi,omega,nu,stream6);
    criticaldt(1,G,theta,phi,omega,nu,streamSD1);               //Stability test for critical timestep
    criticaldt(0.01,G,theta,phi,omega,nu,streamSD2);
    criticaldt(100,G,theta,phi,omega,nu,streamSD3);

    return 0;
}


/*	FUNCTIONS
--------------------------------------------------------------------------------*/

//Energy calculation
double energy(double theta, double omega, double phi, double nu, double R)
{
    return ((omega*omega) + (R*(omega+nu)*(omega+nu)) + ((1+R)*theta*theta) + (R*phi*phi));
}

//Stability testing for timestep dt
void stability(double Ei, double Ef, double dt, ofstream& stream)
{
    if(Ef > 10000000*Ei){ stream << dt << '\t' << 1 << endl; } //Output 1 if unstable
    else{ stream << dt << ' '<< 0 << endl; } //Output 0 if stable
}

//Print t, theta, phi and total energy
void print(double t, double omega, double theta, double phi, double nu, double R, ofstream& stream)
{
    stream << t << ' ' << theta << ' ' << phi << ' ' << energy(theta,omega,phi,nu,R) << endl;
}

//Differentiating Parameters
double ftheta(double omega)
{
    return omega;
}
double fphi(double nu)
{
    return nu;
}
double fomega(double theta, double phi, double omega, double R, double G)
{
    return (-theta*(R+1)+phi*R-omega*G);
}
double fnu(double theta, double phi, double omega, double nu, double R, double G)
{
    return (theta*(R+1)-phi*(R+1)+omega*G*(1-1/R)-nu*G/R);
}

//Print Runge-Kutta solutions
void rungekutta(double R, double G,  int N, double dt, double theta, double phi, double omega, double nu, ofstream& stream)
{
    double t=0, p1=0, p2=0, p3=0, p4=0, k1=0, k2=0, k3=0, k4=0, u1=0, u2=0, u3=0, u4=0, w1=0, w2=0, w3=0, w4=0;
    for(int i=1; i<N; i++)
    {
        print(t,omega,theta,phi,nu,R,stream);

        p1= dt*ftheta(omega);
        k1= dt*fomega(theta,phi,omega,R,G);
        u1= dt*fphi(nu);
        w1= dt*fnu(theta,phi,omega,nu,R,G);

        p2= dt*ftheta(omega+0.5*k1);
        k2= dt*fomega((theta+0.5*p1),(phi+0.5*u1),(omega+0.5*k1),R,G);
        u2= dt*fphi(nu+0.5*w1);
        w2= dt*fnu((theta+0.5*p1),(phi+0.5*u1),(omega+0.5*k1),(nu+0.5*w1),R,G);

        p3= dt*ftheta(omega+0.5*k2);
        k3= dt*fomega((theta+0.5*p2),(phi+0.5*u2),(omega+0.5*k2),R,G);
        u3= dt*fphi(nu+0.5*w2);
        w3= dt*fnu((theta+0.5*p2),(phi+0.5*u2),(omega+0.5*k2),(nu+0.5*w2),R,G);

        p4= dt*ftheta(omega+k3);
        k4= dt*fomega((theta+p3),(phi+u3),(omega+k3),R,G);
        u4= dt*fphi(nu+w3);
        w4= dt*fnu((theta+p3),(phi+u3),(omega+k3),(nu+w3),R,G);

        theta= theta+(p1+2*p2+2*p3+p4)/6;
        omega= omega+(k1+2*k2+2*k3+k4)/6;
        phi= phi+(u1+2*u2+2*u3+u4)/6;
        nu= nu+(w1+2*w2+2*w3+w4)/6;
        t+=dt;
    }
}

//Stability Test for dt
void rungekuttaS(double R, double G,  int N, double dt, double theta, double phi, double omega, double nu, ofstream& stream)
{
    double p1=0, p2=0, p3=0, p4=0, k1=0, k2=0, k3=0, k4=0, u1=0, u2=0, u3=0, u4=0, w1=0, w2=0, w3=0, w4=0;
    double Ei= energy(theta,omega,phi,nu,R), Ef=0;
    for(int i=1; i<N; i++)
    {
        p1= dt*ftheta(omega);
        k1= dt*fomega(theta,phi,omega,R,G);
        u1= dt*fphi(nu);
        w1= dt*fnu(theta,phi,omega,nu,R,G);

        p2= dt*ftheta(omega+0.5*k1);
        k2= dt*fomega((theta+0.5*p1),(phi+0.5*u1),(omega+0.5*k1),R,G);
        u2= dt*fphi(nu+0.5*w1);
        w2= dt*fnu((theta+0.5*p1),(phi+0.5*u1),(omega+0.5*k1),(nu+0.5*w1),R,G);

        p3= dt*ftheta(omega+0.5*k2);
        k3= dt*fomega((theta+0.5*p2),(phi+0.5*u2),(omega+0.5*k2),R,G);
        u3= dt*fphi(nu+0.5*w2);
        w3= dt*fnu((theta+0.5*p2),(phi+0.5*u2),(omega+0.5*k2),(nu+0.5*w2),R,G);

        p4= dt*ftheta(omega+k3);
        k4= dt*fomega((theta+p3),(phi+u3),(omega+k3),R,G);
        u4= dt*fphi(nu+w3);
        w4= dt*fnu((theta+p3),(phi+u3),(omega+k3),(nu+w3),R,G);

        theta= theta+(p1+2*p2+2*p3+p4)/6;
        omega= omega+(k1+2*k2+2*k3+k4)/6;
        phi= phi+(u1+2*u2+2*u3+u4)/6;
        nu= nu+(w1+2*w2+2*w3+w4)/6;

        if (energy(theta,omega,phi,nu,R)>Ef)
        {
            Ef = energy(theta,omega,phi,nu,R);
        }
    }
    stability(Ei,Ef,dt,stream);
}

//Loop function for finding critical timestep
void criticaldt(double R, double G, double theta, double phi, double omega, double nu, ofstream& stream)
{
    int N=0;
    for(double dt=0.001; dt <=5; dt+=0.001)
    {
        N= ceil(150/dt);
        rungekuttaS(R,G,N,dt,theta,phi,omega,nu,stream);
    }
}
