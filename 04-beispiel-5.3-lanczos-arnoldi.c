/*    beispiel-5.3-lanczos-arnoldi.c   02.06.2016            */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<float.h>
#include<limits.h>
#include<string.h>
#include<ctype.h>
#include<sys/time.h>

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif


#define FINT int  /* definiere den Standard F-integer auf der Maschine*/

/* Definiere  Gitterpunkte und Gewichte und  
   Eingabe Parameter als globale Variabeln */

int np=10;                /* Anzahl der Impulsgitterpunkte */
double pmax=4000.0;       /* Maximaler Impuls  beruecksichtigt in MeV */
double V0=-65.1090/197.3; /* Staerke Parameter des Potentials (dimensionslos) */
double mu=0.633*197.3;    /* Reichweite Parameter des Potential in MeV */
double mass=469.34;       /* reduzierte Masse des Systems in MeV */ 
int maxiter=5;            /* maximal Anzahl der Iterationen */
  
double *pmesh,*wmesh;     /* Zeiger auf das Feld mit Gitterpunkten und Gewichten */
double *Vmesh;            /* Zeiger auf ein Feld mit den Potentialmatrixelementen */ 
double *psiwf;            /* Zeiger auf das Feld mit der letzten berechneten Wellenfunktion */

/* dgeev_ als externe Funktion deklarieren 
   Die hier gegebenen Parameter stimmen mit der Definition
   in LAPACK Bibliothek ueberein. 
   Uebergabe der Adresse eines Objekts auch fuer nicht-Felder 
    (in FORTRAN ist das Standard)       
   FORTRAN Name DGEEV    ->    dgeev_                      */

extern void dgeev_(char *jobvl,char *jobvr,FINT *n,double *a,FINT *lda,
                   double *wr,double *wi,double *vl,FINT *ldvl, 
                   double *vr,FINT *ldvr, 
                   double *work,FINT *lwork,FINT *info);


/* Routine, die Gitterpunkte und Gewichte fuer ein Intervall [a,b] festlegt */

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


/* Routine fuer Nullstellensuche mit dem Sekantenverfahren  (siehe Beispiel 3.3-secant.c) */ 
double secant(double x1, double x2, double (*func)(double), int *schritt)
/* x1,x2     Startwerte 
   func      ist die "Referenz" auf eine Funktion mit einem double Parameter, 
             die double zurueckgibt (Referenz = Adresse der Funktion)  
   schritt   ist auch Referenz: Veraenderungen an der Variable wirken sich auf 
                                das aufrufende Programm aus !!! */
{
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



/* =========================== GETUSERPARAM ==================================== */
/*  implemented to set above parameters */

void GetUserParam( int argc, char *argv[] ){

/* Variablen: */
  int i,j;
    char* endptr;
    const char* usage = 
        "mombound [-n <np> -V <V0> -p <pmax> -u <mu> -m <maxiter> -M <mass>]";
    const char* error_message =
        "# FEHLER(GetuserParam): falsche Option: ";

    if (argc>1) { /* falls es ueberhaupt Parameter gibt ... */
        for (i=1; i<argc; i++){
            /* parameter 2 Charakter lang und sollte mit '-' anfaengen ... */
            if ( (strlen(argv[i])==2) && (argv[i][0] == '-') ) { 
                switch (argv[i][1]) { 
                    case 'V': /* falls dies 'V' ist ... */ 
                        /* konvertiere das naechste (++i!) "String"-Argument 
                           nach "double" und weise es der globalen 
                           "double" Variable "V0" zu */
                        V0 = strtod( argv[++i], &endptr); 
                        if ( (!isspace(*endptr) && (*endptr) != 0) ) { 
                            /* sollte mit Leerzeichen abgeschlossen sein */
                            printf(" %s \n %s \n",error_message,usage); 
                            exit(1);
                        }
                        break;
                    case 'M': /* falls dies 'M' ist ... */ 
                        /* konvertiere das naechste (++i!) "String"-Argument 
                           nach "double" und weise es der globalen 
                           "double" Variable "mass" zu */
		        mass = strtod( argv[++i], &endptr); 
                        if ( (!isspace(*endptr) && (*endptr) != 0) ) { 
                            /* sollte mit Leerzeichen abgeschlossen sein */
                            printf(" %s \n %s \n",error_message,usage); 
                            exit(1);
                        }
                        break;

                    case 'p': /* falls dies 'p' ist ... */ 
                        /* konvertiere das naechste (++i!) "String"-Argument 
                           nach "double" und weise es der globalen 
                           "double" Variable "pmax" zu */
                        pmax = strtod( argv[++i], &endptr); 
                        if ( (!isspace(*endptr) && (*endptr) != 0) ) { 
                            /* sollte mit Leerzeichen abgeschlossen sein */
                            printf(" %s \n %s \n",error_message,usage); 
                            exit(1);
                        }
                        break;

                    case 'n':
                        /* konvertiere das naechste (++i!) "String"-Argument 
                           nach "int" und weise es der globalen 
                           "int" Variable "np" zu; 
                           Basis der Konversion ist dezimal (10) */
                        np = strtol( argv[++i], &endptr, 10);
                        if ( (!isspace(*endptr) && (*endptr) != 0) ) {
                            printf(" %s \n %s \n",error_message,usage); 
                            exit(1);
                        }
                        break;
                    case 'm':
                        /* konvertiere das naechste (++i!) "String"-Argument 
                           nach "int" und weise es der globalen 
                           "int" Variable "maxiter" zu; 
                           Basis der Konversion ist dezimal (10) */
                        maxiter = strtol( argv[++i], &endptr, 10);
                        if ( (!isspace(*endptr) && (*endptr) != 0) ) {
                            printf(" %s \n %s \n",error_message,usage); 
                            exit(1);
                        }
                        break;
	
                    default:
                        printf(" %s \n %s \n",error_message,usage); 
                        exit(1);
                }
            } else {
	        printf(" %s \n %s \n",error_message,usage); 
                exit(1);
            } /* end of: if-then-else */
        } /* end-of: for */ 
    } /* end-of: if */

}


/* Funktion, die das Potentialmatrix Element bestimmt 
   input:  p,pp            Impulse in MeV 
   Rueckgabe: Vpot(p,pp)   Potentialmatrixelement in MeV^-2 */

double Vpot(double p, double pp)
  {
    if(p==0.0) p=1e-5;   /* Potential ist nicht definiert fuer p=0 aber */
    if(pp==0.0) pp=1e-5; /* stetig fortgesetzt werden */
    
    return 2.0*V0/(M_PI*4*p*pp)*log(((p+pp)*(p+pp)+mu*mu)/((p-pp)*(p-pp)+mu*mu));
  }  


/* Routine, die die Berechnung der Matrix A_ij vorbereitet
   insbesondere Gitterpunkte und Gewichte festlegt und 
   auch das Potential fuer diese Giterpunkte berechnet.  */ 


void prepmat()   /* Routine benoetigt hier keine Parameter, da diese durch globale Variablen festgelegt werden */
 {
   int i,j;      /* fuer Schleifen ueber Gitterpunkte */

   /* alloziere Speicher fuer Gitterpunkte, Gewichte und das Potential und die Wellenfunktion */

   pmesh=malloc(sizeof(double)*np);
   wmesh=malloc(sizeof(double)*np);
   Vmesh=malloc(sizeof(double)*np*np);
   psiwf=malloc(sizeof(double)*np);
   
   /* lege die Gitterpunkte hier einfach mit bekannter Trapezroutine fest */

   trapez(np,0.0,pmax,pmesh,wmesh);

   /* bestimme Potential an den Stuetzstellen und multipliziere mit Energieunabhaengigen Faktor*/

   for(i=0; i<np; i++)
    for(j=0; j<np; j++)
      {
	Vmesh[i+j*np]=Vpot(pmesh[i],pmesh[j])*pmesh[j]*pmesh[j]*wmesh[j];
      }

 }  
   
/* Einfache Routine zur Addition von Vektoren der Laenge np 
   input:   v1,v2   Zeiger auf Eingabevektoren
   output:  vout    Zeiger auf Ergebnisvektor 
*/

void addvec(double *v1,double *v2, double *vout)
  {
    int i; /* fuer Schleife */

    for(i=0; i<np ; i++)
      {
        vout[i]=v1[i]+v2[i]; /* komponentenweise Addition */
      }    
  }

/* Einfache Routine zur Multiplikation eines Vektors der Laenge np 
   mit einem Skalar
   input:   a       skalarer Faktor 
            vin     Zeiger auf Eingabevektoren
   output:  vout    Zeiger auf Ergebnisvektor 
*/

void multvec(double a,double *vin, double *vout)
  {
    int i; /* fuer Schleife */

    for(i=0; i<np ; i++)
      {
        vout[i]=a*vin[i];      /* komponentenweise Multiplikation */
      }    
  }


/* Einfache Routine zur Berechnung des Skalarproduktes von Vektoren der Laenge np 
   Nutze hier Integration 
   input:       v1,v2   Zeiger auf Eingabevektoren
   Rueckgabe:   v1.v2 
*/

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

/* Routine, die Matrix A(E) auf Vektor anwendet 
   input:   E       Energie in  MeV  
            psiin   Zeiger auf Anfangsvektor(Zustand) (Laenge np in globaler Variable) 
   output:  psiout  Zeiger auf Ergebnisvektor(Zustand) (muss vorher alloziert sein) 
*/

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



/* Routine, die das Lanczos-Arnoldi Verfahren implementiert 
   input:        E Energie in MeV 
   Rueckgabe:    groesster Eigenwert des Problems fuer Energy E 
*/

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

    amat=malloc(sizeof(double)*dim*dim);
    v=malloc(sizeof(double)*dim*np);    
    wvec=malloc(sizeof(double)*np);    
    wtildevec=malloc(sizeof(double)*np);    
    vhilf=malloc(sizeof(double)*np);    
    WR=malloc(sizeof(double)*dim);
    WI=malloc(sizeof(double)*dim);
    work=malloc(sizeof(double)*4*dim);
    c=malloc(sizeof(double)*dim);
    VR=malloc(sizeof(double)*dim*dim);
    VL=malloc(sizeof(double)*dim*dim);
    
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


double func(double E)
  {
    return lambda(E)-1.0;
  }


int main(int argc, char *argv[])
{
  int iterations;
  double ener;
  double timeused;                /* Variabeln zur Zeitmessung */
  struct timeval tvstart,tvend;
 
  /* hier zusaetzlich Routine aus time, um Zeit fuer Loesung zu messen */
  gettimeofday(&tvstart,NULL);

  /* Werte Parameter aus */
  GetUserParam(argc,argv);
  
  /* bereite die Anwendung vor */
  prepmat();

  /* suche Bindungsenergie */ 
  ener=secant(-2.2, -2.0, &func, &iterations);

  printf("Energie: %15.6e MeV in %d Iterationen\n",ener,iterations);

  /* Wie lange hat das in sec gedauert ? */
  gettimeofday(&tvend,NULL);
  timeused=tvend.tv_sec-tvstart.tv_sec;
  timeused=timeused+(tvend.tv_usec-tvstart.tv_usec)*1e-6;  /* Zeitdifferenz  in sec  */

  printf(" \n\n Zeit: %15.6le \n\n",timeused); 

  free(pmesh);
  free(wmesh);
  free(Vmesh);
  free(psiwf);
  
  return 0;
}

/* compilieren und aufrufen mit 
> gcc -llapack -lgfortran beispiel-5.3-lanczos-arnoldi.c -o mombound 
> ./mombound -n 2000 -p 20000 -m 10

 Energie:   -2.208541e+00 MeV in 5 Iterationen
 
 Zeit:    5.393256e+00 
*/
