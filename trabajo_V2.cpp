#include <iostream>
#include <cmath>
#include "Vector.h"
#include "Random64.h"

using namespace std;

const double CHI   =0.1786178958448091;
const double LAMBDA=-0.2123418310626054;
const double XI    =-0.06626458266981849;
const double uno_m2LAMBDA_S2=(1-2*LAMBDA)/2;
const double uno_m2_XIplusCHI=(1-2*(XI+CHI));

const int N=400;
//------------------------- CLASE CUERPO -------------------------------
class Cuerpo{
private:
  double x,Vx,Fx;
  double k,m,R; 
public:
  void Inicie(double x0,double Vx0,double k0,double m0,double R0);
  void CalculeFuerza(double xs);
  void Mueva_r(double dt,double CTE);
  void Mueva_V(double dt,double CTE);
  void Dibujese(void);
  double Getx(void){return x;};
  double GetPx(void){return m*Vx;};
};	
void Cuerpo::Inicie(double x0,double Vx0,double k0,double m0,double R0){
  x=x0; Vx=Vx0; k=k0; m=m0; R=R0;
}
void Cuerpo::CalculeFuerza(double xs){
  Fx=-k*(x-xs);
}
void Cuerpo::Mueva_r(double dt,double CTE){
  x+=Vx*(CTE*dt);
}
void Cuerpo::Mueva_V(double dt,double CTE){
  Vx+=Fx*(CTE*dt/m);
}
void Cuerpo::Dibujese(void){
  cout<<", "<<x<<"+"<<R<<"*cos(t),"<<m*Vx<<"+"<<R<<"*sin(t)";
}
//---------------------- CLASE THERMALBATH ----------------------------
class ThermalBath{
private:
  double KBT, k_resorte, m, Gamma;
  double x0;
public:
  void Inicie(double KBT0, double k_resorte0, double m0, double Gamma0, double x0);
  void Muevase(double dt,Crandom & ran64);
  double Getx(void){return x0;};
};
void ThermalBath::Inicie(double KBT0, double k_resorte0, double m0, double Gamma0, double x00){
  KBT=KBT0; k_resorte=k_resorte0; m=m0; Gamma=Gamma0; x0=x00;
}
void ThermalBath::Muevase(double dt, Crandom & ran64){
  double aux = exp((-2*k_resorte*dt)/(m*Gamma));
  double mu = x0*aux;
  double sigma = sqrt((KBT/k_resorte)*(1-aux));
  x0 = ran64.gauss(mu,sigma);
}

//------------------------- FUNCIONES GLOBALES -------------------------
void InicieAnimacion(double XmaxResorte, double PmaxResorte){
    cout<<"set terminal gif animate"<<endl; 
    cout<<"set output 'Trabajo_V2.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange ["<<-1.8*XmaxResorte<<":"<<1.8*XmaxResorte<<"]"<<endl;
  cout<<"set yrange ["<<-1.8*PmaxResorte<<":"<<1.8*PmaxResorte<<"]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;
}
void InicieCuadro(void){
  cout<<"plot 0,0 ";
}
void TermineCuadro(void){
  cout << endl;
  //	cout << "pause 10" << endl;
}
void GeneraCondicionInicial(double Masa, double Kresorte, double & x0,
			    double & Px0, double XmaxResorte, double PmaxResorte, 
			    double Enmin, double deltaE, Crandom & ran64){
  bool Fuera=true; double E; 
  do{
    x0=(2*ran64.r()-1)*XmaxResorte;
    Px0=(2*ran64.r()-1)*PmaxResorte;
    E=(Px0*Px0/Masa+Kresorte*x0*x0)/2;
    if (E>Enmin && E<(Enmin+deltaE)) Fuera=false;
  }while(Fuera);
}

//------------------------- PROGRAMA PRINCIPAL -------------------------

int main(void){
  //-------------------------Definir e Iniciar Variables
  // EL GENERADOR ALEATORIO
  Crandom ran64(0);
  
  //  LOS RESORTES
  Cuerpo Resorte[N];
  double MResorte=1.0, KResorte=1.0;
  double OmegaResorte=sqrt(MResorte/KResorte);
  double Emin=10, DeltaE=4.0, Xmax=sqrt(2*(Emin+DeltaE)/KResorte), Pmax=sqrt(2*MResorte*(Emin+DeltaE));
  double XResorte,PResorte,VResorte,xs;

   // PARAMETROS DE TIEMPO
  double t,Deltat;
  double T,tmax,tdibujo;
  int ndibujos=100;
  T=2*M_PI/OmegaResorte;  tmax=10*T; Deltat=T/2000;  tdibujo=tmax/ndibujos;
  
  int i;  
  for(i=0; i<N; i++){
    //-----------------(x0,Vx0,k0      ,m0  ,R0);
    GeneraCondicionInicial(MResorte, KResorte, XResorte, PResorte, Xmax, Pmax, Emin, DeltaE, ran64);
    VResorte=PResorte/MResorte;
    Resorte[i].Inicie(XResorte, VResorte, KResorte, MResorte, 0.05);
  }
  
  // EL BANHO TERMICO
  ThermalBath Entorno;
  double KBTBanho=1.0, KBanho=1.0, MBanho=1.0, XBanho0=0.0, GammaBanho=6.0;
  Entorno.Inicie(KBTBanho, KBanho, MBanho, GammaBanho, XBanho0);

  //Termalizando en banho
  for(t=0; t<10.0/GammaBanho; t+=Deltat)
    Entorno.Muevase(Deltat, ran64); 
  
  //CORRER LA SIMULACION
  //Inicie Animación
  InicieAnimacion(Xmax, Pmax);
  
  for(t=0,tdibujo=0; t<tmax; t+=Deltat,tdibujo+=Deltat){
    // Evolucionar el Banho Termico
    Entorno.Muevase(Deltat, ran64);
    xs=Entorno.Getx();
    //Graficar
    if(tdibujo > tmax/ndibujos){
	InicieCuadro();
	for(i=0; i<N; i++){Resorte[i].Dibujese();}
	TermineCuadro();
	tdibujo=0;
      }
    // Evolucionar los resortes
    if(t>tmax/4) //para que inicie quieto
      for(i=0; i<N; i++){
	//Mover el resorte por PEFRL
	Resorte[i].Mueva_r(Deltat,CHI);
	Resorte[i].CalculeFuerza(xs); Resorte[i].Mueva_V(Deltat,uno_m2LAMBDA_S2);
	Resorte[i].Mueva_r(Deltat,XI);
	Resorte[i].CalculeFuerza(xs); Resorte[i].Mueva_V(Deltat,LAMBDA); 
	Resorte[i].Mueva_r(Deltat,uno_m2_XIplusCHI);
	Resorte[i].CalculeFuerza(xs); Resorte[i].Mueva_V(Deltat,LAMBDA); 
	Resorte[i].Mueva_r(Deltat,XI);
	Resorte[i].CalculeFuerza(xs); Resorte[i].Mueva_V(Deltat,uno_m2LAMBDA_S2);
	Resorte[i].Mueva_r(Deltat,CHI);
      }
  }
  
 /* //IMPRIMIR
  	InicieAnimacion(Xmax, Pmax);
  		{
  			InicieCuadro();
  		for(i=0; i<N; i++){Resorte[i].Dibujese();}
  			TermineCuadro();
  			tdibujo=0;
  		}*/
  
  return 0;
}
