// make hausaufgabe4
// ./hausaufgabe4
// Hendrik Falkenberg

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include <gsl/gsl_integration.h>

#include<vector>
#include<iostream>
#include<fstream>

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

#define FINT int  /* definiere den Standard F-integer auf der Maschine*/

int np=400;                /* Anzahl der Impulsgitterpunkte */
double pmax=200.0;       /* Maximaler Impuls  beruecksichtigt in MeV */
double V0=400; /* Staerke Parameter des Potentials (dimensionslos) */
double mu=1;    /* Reichweite Parameter des Potential in MeV */
double mass=1;       /* reduzierte Masse des Systems in MeV */
int maxiter=5;            /* maximal Anzahl der Iterationen */

double *x, *w;
double *pmesh,*wmesh;     /* Zeiger auf das Feld mit Gitterpunkten und Gewichten */
double *Vmesh;            /* Zeiger auf ein Feld mit den Potentialmatrixelementen */
double *psiwf;            /* Zeiger auf das Feld mit der letzten berechneten Wellenfunktion */

void gausslegendre(double a,double b,double *x,double *w,size_t n)
{ gsl_integration_glfixed_table *xwtable;
size_t i;
xwtable=gsl_integration_glfixed_table_alloc(n);
if(xwtable==NULL)
{
printf("Problem with Gauss-Legendre\n");
abort();
}
for(i=0;i<n;i++)
{
gsl_integration_glfixed_point (a, b, i, &x[i], &w[i], xwtable);
}
gsl_integration_glfixed_table_free(xwtable);
}

extern void dgeev_(char *jobvl,char *jobvr,FINT *n,double *a,FINT *lda,
                   double *wr,double *wi,double *vl,FINT *ldvl,
                   double *vr,FINT *ldvr,
                   double *work,FINT *lwork,FINT *info);

double secant(double x1, double x2, double (*func)(double), int *schritt)
{
    /* x1,x2     Startwerte
   func      ist die "Referenz" auf eine Funktion mit einem double Parameter,
             die double zurueckgibt (Referenz = Adresse der Funktion)
   schritt   ist auch Referenz: Veraenderungen an der Variable wirken sich auf
                                das aufrufende Programm aus !!! */
  const double tol=1e-12; /* geforderte Genauigkeit als Konstante */
  double xn;              /* neuer Schaetzwert */

  *schritt=0;  /* noch kein Schritt */

  do
    {         /* naechster Schaetzwert x1,x2 -> xn*/
      xn=x2-func(x2)*(x2-x1)/(func(x2)-func(x1));
      x1=x2;     /* bereite den naechsten Schritt vor:  x2 -> x1 */
      x2=xn;     /*                                     xn -> x2 */

      (*schritt)++;      /* Schritte=Schritte+1 */
    }
  while(fabs(x2-x1)>tol);   /* solange Genauigkeitsziel nicht erreicht */

  return xn;   /* Gebe Nullstelle zurueck */
}

void trapez(int n, double a, double b, double *xp, double *wp)
/* n legt die Anzahl der Stuetzstellen fest.
   a,b sind "normale" double Parameter
   xp und wp sind Zeiger(Pointer) auf ein double Feld,
   es wird die Adresse des Feldes gespeichert !!!  */
{
  int i;
  double h;

  h=(b-a)/(double)(n-1);    /* Berechne Schrittweite */

  for(i=1;i<n-1;i++)
    {
      xp[i]=a+i*h;      /* xp = Anfangsadresse des Feldes */
      wp[i]=h;          /* xp[i] = Nehme die Speicherstelle, */
                        /*      die i Speicherstellen weiter liegt */
    }

  xp[0]=a;          /* Lege Punkte und Gewichte am Rand fest */
  wp[0]=h/2.0;

  xp[n-1]=b;
  wp[n-1]=h/2.0;

}

//void GetUserParam( int argc, char *argv[] ){
//
// /* Variablen: */
//  int i,j;
//    char* endptr;
//    const char* usage =
//        "mombound [-n <np> -V <V0> -p <pmax> -u <mu> -m <maxiter> -M <mass>]";
//    const char* error_message =
//        "# FEHLER(GetuserParam): falsche Option: ";
//
//    if (argc>1) { /* falls es ueberhaupt Parameter gibt ... */
//        for (i=1; i<argc; i++){
//            /* parameter 2 Charakter lang und sollte mit '-' anfaengen ... */
//            if ( (strlen(argv[i])==2) && (argv[i][0] == '-') ) {
//                switch (argv[i][1]) {
//                    case 'V': /* falls dies 'V' ist ... */
//                        /* konvertiere das naechste (++i!) "String"-Argument
//                           nach "double" und weise es der globalen
//                           "double" Variable "V0" zu */
//                        V0 = strtod( argv[++i], &endptr);
//                        if ( (!isspace(*endptr) && (*endptr) != 0) ) {
//                            /* sollte mit Leerzeichen abgeschlossen sein */
//                            printf(" %s \n %s \n",error_message,usage);
//                            exit(1);
//                        }
//                        break;
//                    case 'M': /* falls dies 'M' ist ... */
//                        /* konvertiere das naechste (++i!) "String"-Argument
//                           nach "double" und weise es der globalen
//                           "double" Variable "mass" zu */
//		        mass = strtod( argv[++i], &endptr);
//                        if ( (!isspace(*endptr) && (*endptr) != 0) ) {
//                            /* sollte mit Leerzeichen abgeschlossen sein */
//                            printf(" %s \n %s \n",error_message,usage);
//                            exit(1);
//                        }
//                        break;
//
//                    case 'p': /* falls dies 'p' ist ... */
//                        /* konvertiere das naechste (++i!) "String"-Argument
//                           nach "double" und weise es der globalen
//                           "double" Variable "pmax" zu */
//                        pmax = strtod( argv[++i], &endptr);
//                        if ( (!isspace(*endptr) && (*endptr) != 0) ) {
//                            /* sollte mit Leerzeichen abgeschlossen sein */
//                            printf(" %s \n %s \n",error_message,usage);
//                            exit(1);
//                        }
//                        break;
//
//                    case 'n':
//                        /* konvertiere das naechste (++i!) "String"-Argument
//                           nach "int" und weise es der globalen
//                           "int" Variable "np" zu;
//                           Basis der Konversion ist dezimal (10) */
//                        np = strtol( argv[++i], &endptr, 10);
//                        if ( (!isspace(*endptr) && (*endptr) != 0) ) {
//                            printf(" %s \n %s \n",error_message,usage);
//                            exit(1);
//                        }
//                        break;
//                    case 'm':
//                        /* konvertiere das naechste (++i!) "String"-Argument
//                           nach "int" und weise es der globalen
//                           "int" Variable "maxiter" zu;
//                           Basis der Konversion ist dezimal (10) */
//                        maxiter = strtol( argv[++i], &endptr, 10);
//                        if ( (!isspace(*endptr) && (*endptr) != 0) ) {
//                            printf(" %s \n %s \n",error_message,usage);
//                            exit(1);
//                        }
//                        break;
//
//                    default:
//                        printf(" %s \n %s \n",error_message,usage);
//                        exit(1);
//                }
//            } else {
//	        printf(" %s \n %s \n",error_message,usage);
//                exit(1);
//            } /* end of: if-then-else */
//        } /* end-of: for */
//    } /* end-of: if */
//
//}

double Vpot(double p, double pp)
  {
    if(p==0.0) p=1e-5;   /* Potential ist nicht definiert fuer p=0 aber */
    if(pp==0.0) pp=1e-5; /* stetig fortgesetzt werden */

    return 2.0*V0/(M_PI*4*p*pp)*log(((p+pp)*(p+pp)+mu*mu)/((p-pp)*(p-pp)+mu*mu));
  }

void prepmat()   /* Routine benoetigt hier keine Parameter, da diese durch globale Variablen festgelegt werden */
 {
   int i,j;      /* fuer Schleifen ueber Gitterpunkte */

   /* alloziere Speicher fuer Gitterpunkte, Gewichte und das Potential und die Wellenfunktion */

   double *pmesh=(double*) malloc(sizeof(double)*np);
   double *wmesh=(double*) malloc(sizeof(double)*np);
   double *Vmesh=(double*) malloc(sizeof(double)*np*np);
   double *psiwf=(double*) malloc(sizeof(double)*np);

   /* lege die Gitterpunkte hier einfach mit bekannter Trapezroutine fest */

   trapez(np,0.0,pmax,pmesh,wmesh);

   /* bestimme Potential an den Stuetzstellen und multipliziere mit Energieunabhaengigen Faktor*/

   for(i=0; i<np; i++)
    for(j=0; j<np; j++)
      {
	Vmesh[i+j*np]=Vpot(pmesh[i],pmesh[j])*pmesh[j]*pmesh[j]*wmesh[j];
      }

 }

 void addvec(double *v1,double *v2, double *vout)
  {
    int i; /* fuer Schleife */

    for(i=0; i<np ; i++)
      {
        vout[i]=v1[i]+v2[i]; /* komponentenweise Addition */
      }
  }

  void multvec(double a,double *vin, double *vout)
  {
    int i; /* fuer Schleife */

    for(i=0; i<np ; i++)
      {
        vout[i]=a*vin[i];      /* komponentenweise Multiplikation */
      }
  }

  double skalp(double *v1,double *v2)
  {
    int i;       /* fuer Schleife */
    double sum;  /* Zwischenergebniss Summation */

    sum=0.0;
    for(i=0; i<np ; i++)
      {
        sum+=v1[i]*v2[i]*pmesh[i]*pmesh[i]*wmesh[i]; /* summation ueber alle Stuetzstellen */
      }
    return sum;
  }

  void applymat(double *psiin,double *psiout,double E)
  {
    int i,j;      /* fuer Schleifen ueber Gitterpunkte */
    double Efakt; /* Variabel fuer energieabhaengigen Faktor */

    for(i=0; i<np; i++)
      {
       psiout[i]=0.0;
       Efakt=2.0*mass/(2*mass*E-pmesh[i]*pmesh[i]);

       for(j=0; j<np; j++)
        {
	  psiout[i]+=Efakt*Vmesh[i+np*j]*psiin[j];

        }
      }

  }

  double lambda(double E)
  {
    int i,j,k;
    double *amat;       /* fuer das Feld A_ij */
    double norm,prod;   /* fuer skalare Zwichenergebnisse */
    double *v;          /* Feld mit Basis-Vektoren */
    double *wvec;       /* Vektor w_{j+1} */
    double *wtildevec;  /* Vektor wtilde_{j+1} */
    double *vhilf;      /* Hilfsvektor fuer Zwischenergebnis */
    double *VL,*VR;     /* Speicher Problems */
    double *c;          /* Eigenvektor des reduzierten Problems */
    double *WR,*WI;     /* Felder fuer Real- und Imaginaerteil der Eigenwerte */
    double *work;       /* Hilfsfeld, kann uns egal sein was dort gespeichert wird */
    double maxlambda;   /* fuer den maximalen Eigenwert */

    /* INTEGER Variabel in FORTRAN = long = 4 Bytes */
    FINT info;       /* info sollte Null sein beim Verlassen der Routine, sonst Fehler */
    FINT lda=maxiter+1,ldvl=maxiter+1,ldvr=maxiter+1;   /* Dimension der Felder */
    FINT lwork=4*(maxiter+1);   /* Dimension des work-Felder (sollte mind. 4*dim sein)*/
    FINT dim=maxiter+1;         /* Dimension des reduzierten Problems */
    char jobvl='N',jobvr='V';   /* char Variabeln bestimmen, ob VR und VL berechnet werden
				 sollen -> hier berechnen, 'N' heisst nicht berechnen */
    /* alloziere Speicher fuer Matrizen */

    amat=(double*)malloc(sizeof(double)*dim*dim);
    v=(double*)malloc(sizeof(double)*dim*np);
    wvec=(double*)malloc(sizeof(double)*np);
    wtildevec=(double*)malloc(sizeof(double)*np);
    vhilf=(double*)malloc(sizeof(double)*np);
    WR=(double*)malloc(sizeof(double)*dim);
    WI=(double*)malloc(sizeof(double)*dim);
    work=(double*)malloc(sizeof(double)*4*dim);
    c=(double*)malloc(sizeof(double)*dim);
    VR=(double*)malloc(sizeof(double)*dim*dim);
    VL=(double*)malloc(sizeof(double)*dim*dim);

    /* belege amat mit Null */


    for(i=0; i<(maxiter+1)*(maxiter+1); i++)
      {
	amat[i]=0.0;
      }

    /* setze Startvektor auf konstant 1 und normiere */

    for(i=0; i<np; i++)
      {
	v[i+0*np]=1.0;
      }

    norm=sqrt(skalp(&v[0*np],&v[0*np])); /* Bestimme Norm^-1 */
    multvec(1.0/norm,&v[0*np],&v[0*np]);         /* multipliziere mit norm^-1 */

    /* Iteriere hier immer bis zur maximalen Iteration
       moegliche Verbesserung: teste auf Konvergenz und iteriere weniger falls moeglich */

    for(k=0; k<maxiter; k++)
      {
	applymat(&v[k*np],wvec,E);              /* w=A*v[k] */

	/* Gram-Schmidt Orthogonalisierung */
	multvec(1.0,wvec,wtildevec);            /* wtilde=wvec */


        for(j=0; j<=k; j++)
	  {
	    prod=skalp(&v[np*j],wvec);          /* prod = (v[j],w) */

	    amat[j+dim*k]=prod;                  /* Zwischenergebnis: A[j,k] = <v[j] | A | v[k] > */
	    multvec(-prod,&v[np*j],vhilf);      /* wtilde=wtilde-v[j]*(v[j],w) */
	    addvec(wtildevec,vhilf,wtildevec);
	  }

	amat[k+1+dim*k]=sqrt(skalp(wtildevec,wtildevec)); /* wegen skalierung die folgt, A[k+1,k]= sqrt(<wtilde,wtilde>) */
        norm=1.0/amat[k+1+dim*k];                     /* Bestimme Norm^-1 */
        multvec(norm,wtildevec,&v[(k+1)*np]);            /* v[k+1] = wtilde / |wtilde| */
      }

    /* hier ist Matrix A bestimmt */
    /* Eigenwerte koenne beispielsweise mit DGEEV bestimmt werden */


    dim--;  /* Achtung letzte Iteration erzeugt nur A[maxiter,maxiter-1] aber kein Element der Spalte A[*,maxiter] */
            /* deswegen dimension niedriger waehlen als Speicher (LDxx != N !!!!) */
    dgeev_(&jobvl,&jobvr,&dim,amat,&lda,WR,WI,VL,&ldvl,VR,&ldvr,work,&lwork,&info);


    maxlambda=0.0;	   /* sehr einfach Maximumsuche, beachtet nur reelle Eigenwerte */
    for(i=0; i<dim; i++)
      {
	if(abs(WI[i])<1e-6) maxlambda=fmax(maxlambda,WR[i]);
      }


    /* Test ob Erfolgreich */
    if(info != 0)
      {
	printf("Info: %d \n",info);
	abort();
      }



    /* Speicher wieder freigeben */

    free(amat);
    free(v);
    free(wvec);
    free(wtildevec);
    free(vhilf);
    free(WR);
    free(WI);
    free(work);
    free(c);
    free(VL);
    free(VR);

    return maxlambda;
  }

double f(double x)
{
return 0;
}

double v_mu_l0 (double &p, double &p_prime)
 {

    double integral_part, result;
    integral_part = gausslegendre(0,10,x,w,100);
    result = 2/M_PI*v0*integral_part;
    return result;
 }

//berechnet die Funktion analytisch
double v_mu_l0_analytic (double &p, double &p_prime, double &mu) {
return 0;
}

void perform_3() {
    //Parameter: a=0, b=200, n=400;
    double a=0;
    double b=200;
    double n=400;
    double step=(b-a)/n;
    double result;

    for(int p=a; p<b; p+=step) {
        for (p_prime =a; p_prime<b; p+=step) {
            v_mu_l0(p, p_prime);
        }

    }


    maxreldiff;

}

void perform_4() {

//    int iterations;
//  double ener;
//  double timeused;                /* Variabeln zur Zeitmessung */
//  struct timeval tvstart,tvend;
//
//  /* hier zusaetzlich Routine aus time, um Zeit fuer Loesung zu messen */
//  gettimeofday(&tvstart,NULL);
//
//  /* Werte Parameter aus */
//  GetUserParam(argc,argv);
//
//  /* bereite die Anwendung vor */
//  prepmat();
//
//  /* suche Bindungsenergie */
//  ener=secant(-2.2, -2.0, &func, &iterations);
//
//  printf("Energie: %15.6e MeV in %d Iterationen\n",ener,iterations);
//
//  /* Wie lange hat das in sec gedauert ? */
//  gettimeofday(&tvend,NULL);
//  timeused=tvend.tv_sec-tvstart.tv_sec;
//  timeused=timeused+(tvend.tv_usec-tvstart.tv_usec)*1e-6;  /* Zeitdifferenz  in sec  */
//
//  printf(" \n\n Zeit: %15.6le \n\n",timeused);
}

void perform_5() {
}

void perform_6() {
}

void perform_7() {
}

void perform_8() {
}

int main() {

    perform_3;
    perform_4;
    perform_5;
    perform_6;
    perform_7;
    perform_8;

    return 0
}
