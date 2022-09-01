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
  void Dibujese(double x0,double Vx0, double w, double t);
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
void Cuerpo::Dibujese(double x0,double Vx0, double w, double t){
	double B = Vx0/w;
  cout<<(x0)*cos(t*w)+(B)*sin(t*w) <<" "<<(((-1*x0))*sin(t*w)+(B)*cos(t*w))*w*m << endl;
}


//------------------------- FUNCIONES GLOBALES -------------------------
void InicieAnimacion(double XmaxResorte, double PmaxResorte){
  cout<<"set terminal pngcairo"<<endl; 
  cout<<"set output 'Resorte.png'"<<endl;
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


//------------------------- PROGRAMA PRINCIPAL -------------------------

int main(void)
{
  double MResorte=1.0, KResorte=1.0;
  double OmegaResorte=sqrt(MResorte/KResorte);
  double Emin=10, DeltaE=4.0, Xmax=sqrt(2*(Emin+DeltaE)/KResorte), Pmax=sqrt(2*MResorte*(Emin+DeltaE));
  double XResorte=1.5,PResorte=1.0,VResorte=PResorte/MResorte;
	double xs=0;
  Cuerpo Resorte;

  // PARAMETROS DE TIEMPO
  double t,Deltat;
  double T,tmax,tdibujo;
  int ndibujos=10000;
  T=2*M_PI/OmegaResorte;  tmax=10*T; Deltat=T/2000;  tdibujo=tmax/ndibujos;

  Resorte.Inicie(XResorte, VResorte, KResorte, MResorte, 0.05);
  //InicieAnimacion(XResorte, PResorte);

  for(t=0,tdibujo=0; t<T; t+=Deltat,tdibujo+=Deltat){
    //Graficar
    if(tdibujo > tmax/ndibujos){
	//InicieCuadro();
	Resorte.Dibujese(XResorte, VResorte, OmegaResorte, t);
	//TermineCuadro();
	tdibujo=0;
      }
    //Mover el resorte por PEFRL
	Resorte.Mueva_r(Deltat,CHI);
	Resorte.CalculeFuerza(xs); Resorte.Mueva_V(Deltat,uno_m2LAMBDA_S2);
	Resorte.Mueva_r(Deltat,XI);
	Resorte.CalculeFuerza(xs); Resorte.Mueva_V(Deltat,LAMBDA); 
	Resorte.Mueva_r(Deltat,uno_m2_XIplusCHI);
	Resorte.CalculeFuerza(xs); Resorte.Mueva_V(Deltat,LAMBDA); 
	Resorte.Mueva_r(Deltat,XI);
	Resorte.CalculeFuerza(xs); Resorte.Mueva_V(Deltat,uno_m2LAMBDA_S2);
	Resorte.Mueva_r(Deltat,CHI);
  }

  return 0;
}

