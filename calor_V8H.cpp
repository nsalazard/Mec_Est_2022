#include <iostream>
#include <cmath>
#include <algorithm>
#include <fstream>
#include "Vector.h"
#include "Random64.h"
#include <set>

using namespace std;

const double CHI   =0.1786178958448091;
const double LAMBDA=-0.2123418310626054;
const double XI    =-0.06626458266981849;
const double uno_m2LAMBDA_S2=(1-2*LAMBDA)/2;
const double uno_m2_XIplusCHI=(1-2*(XI+CHI));
const double HBAR =1.0;

const int Nmax=810;
//------------------------- CLASE CUERPO -------------------------------
class Cuerpo{
private:
  double x,Vx,Fx;
  double k,m,R,GammaResorte; 
public:
  void Inicie(double x0,double Vx0,double k0,double m0,double R0,double GammaResorte0);
  void CalculeFuerza(double xs);
  void Mueva_r(double dt,double CTE);
  void Mueva_V(double dt,double CTE);
  void Dibujese(void);
  double Getx(void){return x;};
  double GetPx(void){return m*Vx;};
};	
void Cuerpo::Inicie(double x0,double Vx0,double k0,double m0,double R0,double GammaResorte0){
  x=x0; Vx=Vx0; k=k0; m=m0; R=R0; GammaResorte=GammaResorte0;
}
void Cuerpo::CalculeFuerza(double xs){
  Fx=-k*(x-xs)-m*GammaResorte*Vx;
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
//---------------------- CLASE RANDOM THERMALBATH ----------------------------
class ThermalBath{
private:
  double Sigma;
public:
  void Inicie(double KBT0, double k_resorte0);
  double Getx(Crandom & ran64);
};
void ThermalBath::Inicie(double KBT0, double k_resorte0){
  Sigma=sqrt(KBT0/k_resorte0);
}
double ThermalBath::Getx(Crandom & ran64){
  return ran64.gauss(0,Sigma);
}

//------------------------- FUNCIONES GLOBALES -------------------------
void InicieAnimacion(double XmaxResorte, double PmaxResorte){
    cout<<"set terminal gif animate"<<endl; 
    cout<<"set output 'calor_60000.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange ["<<-5.0*XmaxResorte<<":"<<5.0*XmaxResorte<<"]"<<endl;
  cout<<"set yrange ["<<-5.0*PmaxResorte<<":"<<5.0*PmaxResorte<<"]"<<endl;
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
  // EL NUMERO INICIAL DE OSCILADORES
  // EL GENERADOR ALEATORIO
  Crandom ran64(0);
  
  //  LOS RESORTES
  double KBT = 1.0 ;
  Cuerpo Resorte[Nmax];
  double MResorte=1.0, KResorte=1.0, GammaResorte=4.0/30000; // GAMMARESORTE!!!!
  double OmegaResorte=sqrt(MResorte/KResorte);
  double Emin=10, DeltaE=4.0,
    Xmax=sqrt(2*(Emin+DeltaE)/KResorte), Pmax=sqrt(2*MResorte*(Emin+DeltaE));
  double XResorte,PResorte,VResorte,xs;
  int i;
  int N=30;
	set<double> histop;
	set<double> histox;
  for(i=0; i<N; i++){
    //-----------------(x0,Vx0,k0      ,m0  ,R0);
    GeneraCondicionInicial(MResorte, KResorte, XResorte, PResorte,
			   Xmax, Pmax, Emin, DeltaE, ran64);
    VResorte=PResorte/MResorte;
    Resorte[i].Inicie(XResorte, VResorte, KResorte, MResorte, 0.05, GammaResorte);
  }
  
  // EL BANHO TERMICO
  ThermalBath Entorno;
  double KBTBanho=1.0, KBanho=1.0;
  Entorno.Inicie(KBTBanho, KBanho);
  
  // PARAMETROS DE TIEMPO
  double t,Deltat;
  double T,tmax,tdibujo;
  int c_multiplicar,N_multiplicar=3, cM,NM=3;
  int ndibujos=100;
  T=2*M_PI/OmegaResorte;  tmax=1000*T; Deltat=T/1000;  tdibujo=tmax/ndibujos;
  /*//CORRER LA SIMULACION
  //Inicie Animaci??n
  InicieAnimacion(Xmax, Pmax);
    //Graficar est??tico para ver la condici??n inicial microcan??nica est??tica
  for(t=0,tdibujo=0; t<tmax; t+=Deltat,tdibujo+=Deltat){
    //Graficar
    if(tdibujo > tmax/ndibujos){
      InicieCuadro();
      for(i=0; i<N; i++){Resorte[i].Dibujese();}
      TermineCuadro();
      tdibujo=0;
    }
  }
  //AHORA S??, CORRE LA SIMULACI??N 
  //Multiplicar los resortes */
  for(c_multiplicar=0;c_multiplicar<N_multiplicar;c_multiplicar++){
    for(cM=1;cM<NM;cM++)
      for(i=0; i<N; i++)
	Resorte[i+cM*N]=Resorte[i];
    N=NM*N;
    
    //Correr tmax con ese n??mero de resortes
    for(t=0,tdibujo=0; t<tmax; t+=Deltat,tdibujo+=Deltat){
      
			/*//Graficar
      if(tdibujo > tmax/ndibujos){
	InicieCuadro();
	for(i=0; i<N; i++){Resorte[i].Dibujese();}
	TermineCuadro();
	tdibujo=0;
      }*/
      
      // Evolucionar los resortes
      if(t>0){ //para que inicie quieto
	for(i=0; i<N; i++){
	  // Para cada resorte calcular el valor de xs impuesto por el banho termico
	  xs=Entorno.Getx(ran64);
	  //Con este valor, mover el resorte por PEFRL
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
		if (t > (tmax -(20*Deltat)))	{
  	for(i=0; i<Nmax; i++){
			histox.insert(Resorte[i].Getx());
			histop.insert(Resorte[i].GetPx());
		}
		}
	
    }
  }
  //IMPRIMIR
    //Calcular la desviaci??n est??ndar en x
  double xprom,x2prom,sigma_x;
  //Calculo xprom
  for(xprom=0,i=0; i<N; i++)
    xprom+=Resorte[i].Getx();
  xprom/=N;
  //Calculo x2prom
  for(x2prom=0,i=0; i<N; i++)
    x2prom+=pow(Resorte[i].Getx(),2.0);
  x2prom/=N;
  //Calculo sigma_x
  sigma_x=sqrt(x2prom-xprom*xprom);
  //Imprimir sigma_x
  //cout<<"sigma_x(computacional)="<<sigma_x<<endl;
  
  //Calcular la desviaci??n est??ndar en xp
  double pprom,p2prom,sigma_p;
  //Calculo xprom
  for(pprom=0,i=0; i<N; i++)
    pprom+=Resorte[i].GetPx();
  pprom/=N;
  //Calculo x2prom
  for(p2prom=0,i=0; i<N; i++)
    p2prom+=pow(Resorte[i].GetPx(),2.0);
  p2prom/=N;
  //Calculo sigma_p
  sigma_p=sqrt(p2prom-pprom*pprom);
  //Imprimir sigma_p
  //cout<<"sigma_p(computacional)="<<sigma_p<<endl;


	ofstream HistoX;
	ofstream HistoP;
	ofstream Prob;
  	HistoX.open ("HistoX.txt");
	HistoP.open ("HistoP.txt");
	Prob.open ("Prob.txt");
	
	double sigmax = 3.37998/10;
	double sigmap =3.55288/10;

	double deltaX = (*(histox.end()) - *(histox.begin()))/sigmax;
	double xmin, xmax, pmin,pmax,ProbX, ProbPx;
	double countx, countp;
	double maxCx= 0, maxCp=0;

	for(int ii=0; ii <= deltaX+1;ii++){
		countx = 0; countp = 0;
		xmin = *(histox.begin())+(sigmax*ii); xmax = *(histox.begin())+(sigmax*(ii+1));
		pmin = *(histop.begin())+(sigmap*ii); pmax = *(histop.begin())+(sigmap*(ii+1));

		
		for(auto ix=histox.begin(),ip=histop.begin(); ix !=histox.end() ; ix++, ip++){
			
			ProbX = sqrt(2*M_PI*MResorte*KBT)*exp((-1)*pow( *ix,2)/(2*KBT/ KResorte ))/(HBAR*(1/KBT)*OmegaResorte);
			ProbPx =sqrt(2*M_PI*KBT/KResorte)*exp((-1)*pow( *ip,2)/(2*MResorte*KBT))/(HBAR*(1/KBT)*OmegaResorte) ;
			Prob << *ix  << "\t" << ProbX/2.5 << "\t" << *ip  << "\t" << ProbPx/2.5 << endl;
			
			if(*ix >= xmin && *ix < xmax){countx += 1;}
			if(*ip >= pmin && *ip < pmax){countp += 1;}
				}

		
		HistoX << ((xmin+xmax)/2) << "\t" << (countx/(39.0)) << endl;
		HistoP << ((pmin+pmax)/2) << "\t" << (countp/(43.0)) << endl;

		if (countx > maxCx){
			maxCx = countx;
		}
		if (countp > maxCp){
			maxCp = countp;
		}
	}
	HistoX.close();
	HistoP.close();
	Prob.close();

	// cout<<"Maximo valor de X"<< maxCx<<endl;
	// cout<<"Maximo valor de P"<<maxCp<<endl;

	// Histograma X
	cout<<"set term pdf"<<endl; 
  cout<<"set out 'HistogramaPosicion.pdf'"<<endl;
	cout<<"set title 'Histograma de la Posicion'"<<endl;
  cout<<"set ylabel 'Numero de particulas'"<<endl;
	cout<<"set xlabel 'Posicion x'"<<endl;
	cout<<"set autoscale"<<endl;
	cout<<"set key"<<endl;
	cout<<"set font ',7'"<<endl;
	cout<<"plot 'HistoX.txt' u 1:2 with impulses t 'Data' ,'Prob.txt' u 1:2 w p t 'Distribucion marginal P(x)'"<<endl;

	// Histograma Y
	cout<<"set term pdf"<<endl; 
  cout<<"set out 'HistogramaMomento.pdf'"<<endl;
	cout<<"set title 'Histograma del Momento'"<<endl;
  cout<<"set ylabel 'Numero de particulas'"<<endl;
	cout<<"set xlabel 'Momento Px'"<<endl;
	cout<<"set autoscale"<<endl;
	cout<<"set key"<<endl;
	cout<<"set font ',7'"<<endl;
	cout<<"plot 'HistoP.txt' u 1:2 with impulses t 'Data','Prob.txt' u 3:4 w p t 'Distribucion marginal P(p)'"<<endl;
	
  return 0;
}
