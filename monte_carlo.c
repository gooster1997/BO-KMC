
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <unistd.h>
#include <time.h>
#include <Python.h>

static void calc_k(double params[], double k[], int VNa){ // calculate Monte Carlo rate constants

    k[0]=params[0];
    k[1]=params[1]/(VNa);
    k[2]=params[2]/(VNa);
    k[3]=params[3]/(VNa);
    k[4]=params[4]/(VNa);
    k[5]=params[5]/(VNa);
    k[6]=params[6]/(VNa);
    k[7]=params[7]/(VNa);
    k[8]=2.0*params[8]/(VNa);
    k[9]=2.0*params[9]/(VNa);
    k[10]=params[10]/(VNa);
    k[11]=params[11]/(VNa);
    k[12]=params[12]/(VNa);
    k[13]=params[13]/(VNa);

    return;
    }

static void calc_r(int Xi1, int Xi2, int Xa, int Xb, long Xatol, long Xbtol, double k[], double r[]){ // calculate reaction rates
    
    r[0]=Xi1*k[0];
    r[1]=Xi2*k[1]*Xa;
    r[2]=Xi2*k[2]*Xb;
    r[3]=Xa*Xatol*k[3];
    r[4]=Xb*Xatol*k[4];
    r[5]=Xa*Xbtol*k[5];
    r[6]=Xb*Xbtol*k[6];
    r[7]=Xatol*Xbtol*k[7];
    r[8]=Xatol*(Xatol-1)*k[8]/2;
    r[9]=Xbtol*(Xbtol-1)*k[9]/2;
    r[10]=Xa*Xatol*k[10];
    r[11]=Xb*Xatol*k[11];
    r[12]=Xa*Xbtol*k[12];
    r[13]=Xb*Xbtol*k[13];
    r[14]=0;
    
    for (int j=0; j<=13; j++) {r[14]=r[14]+r[j];}
    
    return;
    }

static int gene_mu(double r[]){ // generate random number to select reaction
    double r1;
    r1=rand()/(RAND_MAX+1.0);
    double sum=0;
    int mu=0;
    double temp;
    temp=r1*r[14];
    while(sum<temp){
        sum=sum+r[mu];
        mu++;
        }
    return mu;
    }

static int chain_judge(long Xtol){ // generate random number to select chain
    double r=0;
    r=rand()/(RAND_MAX+1.0);
    int num;
    num=Xtol*r;

    return num;
    }

static void homopro(long *X, long Xtol){ // homopropagation
    int num;
    num=chain_judge(Xtol);
    (*(X+num))++;

    return;
    }

static int chain_delete(long *X, int num, long Xtol){ // delete chain for crosspropagation
    int n;
    n=*(X+num);
    for (int i=num;i<=Xtol-2;i++) {*(X+i)=*(X+i+1);}
    return n;
    }

static void crosspro(long *Aa, long *Ab, long *Ba, long *Bb, long Xatol, long Xbtol){ // crosspropagation
    int num, a, b;
    num=chain_judge(Xatol);
    a=chain_delete(Aa, num, Xatol);
    b=chain_delete(Ab, num, Xatol);
    *(Ba+Xbtol)=a;
    *(Bb+Xbtol)=b+1;

    return;
    }

static void homoter(long *a, long *b, long *Pa, long *Pb, long Xtol, int Xptol){ // homotermination
    int num1, num2, a1, a2, b1, b2;
    num1=chain_judge(Xtol);
    a1=chain_delete(a, num1, Xtol);
    b1=chain_delete(b, num1, Xtol);
    num2=chain_judge(Xtol-1);
    a2=chain_delete(a, num2, Xtol-1);
    b2=chain_delete(b, num2, Xtol-1);
    *(Pa+Xptol)=a1+a2;
    *(Pb+Xptol)=b1+b2;

    return;
    }

static void crosster(long *Aa, long *Ab, long *Ba, long *Bb, long *Pa, long *Pb, long Xatol, long Xbtol, int Xptol){ // crosstermination
    int num1, num2, a1, a2, b1, b2;
    num1=chain_judge(Xatol);
    num2=chain_judge(Xbtol);
    a1=chain_delete(Aa, num1, Xatol);
    b1=chain_delete(Ab, num1, Xatol);
    a2=chain_delete(Ba, num2, Xbtol);
    b2=chain_delete(Bb, num2, Xbtol);
    *(Pa+Xptol)=a1+a2;
    *(Pb+Xptol)=b1+b2;

    return;
    }

static void trn(long *A, long *B, long *Pa, long *Pb, long Xtol, int Xptol){ // chain transfer
    int num, a, b;
    num=chain_judge(Xtol);
    a=chain_delete(A, num, Xtol);
    b=chain_delete(B, num, Xtol);
    *(Pa+Xptol)=a;
    *(Pb+Xptol)=b;

    return;
    }

static void reaction_execute(int mu, int *Xi1, int *Xi2, long *Xa, long *Xb, long *Xatol, long *Xbtol, int *Xptol, long *Aa, long *Ab, long *Ba, long *Bb, long *Pa, long *Pb){

    switch(mu){
        case 1:
            (*Xi1)--;
            (*Xi2)+=2;
            break;
        case 2:
            (*Xi2)--;
            (*Xa)--;
            *(Aa+*Xatol)=1;
            *(Ab+*Xatol)=0;
            (*Xatol)++;
            break;
        case 3:
            (*Xi2)--;
            (*Xb)--;
            *(Bb+*Xbtol)=1;
            *(Ba+*Xbtol)=0;
            (*Xbtol)++;
            break;
        case 4:
            homopro(Aa, *Xatol);
            (*Xa)--;
            break;
        case 5:
            crosspro(Aa, Ab, Ba, Bb, *Xatol, *Xbtol);
            (*Xb)--;
            (*Xatol)--;
            (*Xbtol)++;
            break;
        case 6:
            crosspro(Bb, Ba, Ab, Aa, *Xbtol, *Xatol);
            (*Xa)--;
            (*Xbtol)--;
            (*Xatol)++;
            break;
        case 7:
            homopro(Bb, *Xbtol);
            (*Xb)--;
            break;
        case 8:
            crosster(Aa, Ab, Ba, Bb, Pa, Pb, *Xatol, *Xbtol, *Xptol);
            (*Xatol)--;
            (*Xbtol)--;
            (*Xptol)++;
            break;
        case 9:
            homoter(Aa, Ab, Pa, Pb, *Xatol, *Xptol);
            (*Xatol)-=2;
            (*Xptol)++;
            break;
        case 10:
            homoter(Ba, Bb, Pa, Pb, *Xbtol, *Xptol);
            (*Xbtol)-=2;
            (*Xptol)++;
            break;
        case 11:
            trn(Aa, Ab, Pa, Pb, *Xatol, *Xptol);
            (*Xa)--;
            (*Xi2)++;
            (*Xatol)--;
            (*Xptol)++;
            break;
        case 12:
            trn(Aa, Ab, Pa, Pb, *Xatol, *Xptol);
            (*Xb)--;
            (*Xi2)++;
            (*Xatol)--;
            (*Xptol)++;
            break;
        case 13:
            trn(Ba, Bb, Pa, Pb, *Xbtol, *Xptol);
            (*Xa)--;
            (*Xi2)++;
            (*Xbtol)--;
            (*Xptol)++;
            break;
        case 14:
            trn(Ba, Bb, Pa, Pb, *Xbtol, *Xptol);
            (*Xb)--;
            (*Xi2)++;
            (*Xbtol)--;
            (*Xptol)++;
            break;
        }
    }

static double calc_Mn(long *Mtemp, long *Aa, long *Ab, long *Ba, long *Bb, long *Pa, long *Pb, long Xatol, long Xbtol, int Xptol){
    double Mn;

    for (int i=0;i<Xatol;i++){*Mtemp=*Mtemp+120*(*(Aa+i))+200*(*(Ab+i));}
    for (int i=0;i<Xbtol;i++){*Mtemp=*Mtemp+120*(*(Ba+i))+200*(*(Bb+i));}
    for (int i=0;i<Xptol;i++){*Mtemp=*Mtemp+120*(*(Pa+i))+200*(*(Pb+i));}

    Mn=1.0*(*Mtemp)/(Xatol+Xbtol+Xptol);
    return Mn;
    }

static double calc_MWD(long Mtemp, double Mn, long *Aa, long *Ab, long *Ba, long *Bb, long *Pa, long *Pb, long Xatol, long Xbtol, int Xptol){
    long long M2;
    double MWD, Mw;
    M2=0;

    for (int i=0;i<Xatol;i++){M2=M2+(120*(*(Aa+i))+200*(*(Ab+i)))*(120*(*(Aa+i))+200*(*(Ab+i)));}
    for (int i=0;i<Xbtol;i++){M2=M2+(120*(*(Ba+i))+200*(*(Bb+i)))*(120*(*(Ba+i))+200*(*(Bb+i)));}
    for (int i=0;i<Xptol;i++){M2=M2+(120*(*(Pa+i))+200*(*(Pb+i)))*(120*(*(Pa+i))+200*(*(Pb+i)));}

    Mw=1.0*M2/Mtemp;
    MWD=1.0*Mw/Mn;
    return MWD;
    }

static void KMC(double params[], double *Mn, double *MWD, double *conv){ // main function
	// Assign parameters
	double k[14] = {0}; 
	double conc, ratio; 
	int VNa = 50000000;
	calc_k(params, k, VNa); // calculate Monte Carlo rate constants
	ratio = params[14]; //ratio of initiator and monomer
	conc = params[15]; // monomer concentration

    // Initialization	
	int Xi1, Xi2, Xptol, Xi1thre;
	long Xatol, Xbtol, Xa, Xb;
	Xa = Xb = VNa * conc; // number of comonomers
	Xi1 = Xa * ratio; // number of initiator
	Xi1thre = 0.5 * Xi1; // threshold number of initiator for triggering termination of simulation
	Xi2 = Xptol = 0; // Xi2: number of primary radical, Xptol: number of dead polymer
	Xatol = Xbtol = 0; // number of living radicals

    long *Aa, *Ab, *Ba, *Bb, *Pa, *Pb;
    Aa = (long *)malloc(10000000 * sizeof(long)); // number of monomer a in radicals ended with monomer a
    Ab = (long *)malloc(10000000 * sizeof(long)); // number of monomer b in radicals ended with monomer a
    Ba = (long *)malloc(10000000 * sizeof(long)); // number of monomer a in radicals ended with monomer b
    Bb = (long *)malloc(10000000 * sizeof(long)); // number of monomer b in radicals ended with monomer b
    Pa = (long *)malloc(100000000 * sizeof(long)); // number of monomer a in dead polymers
    Pb = (long *)malloc(100000000 * sizeof(long)); // number of monomer b in daed polymers
	double r[15]={0}; // reaction rates
  
    int mu; // index of selected reaction 
    long cycle=1;
    long Mtemp;

    srand(2021);
	time_t start;
	start = time(NULL);

	// Start iteration
    while((Xa+Xb)>=0 && Xi1>Xi1thre){
        calc_r(Xi1, Xi2, Xa, Xb, Xatol, Xbtol, k, r); // calculate reaction rates
        mu = gene_mu(r); // randomly select reaction to execute
        reaction_execute(mu, &Xi1, &Xi2, &Xa, &Xb, &Xatol, &Xbtol, &Xptol, Aa, Ab, Ba, Bb, Pa, Pb);
        
        if((time(NULL)-start)>900) {break;}
        cycle++;
        }
	
	// Calculate final results
    Mtemp = 0;
    *Mn = calc_Mn(&Mtemp, Aa, Ab, Ba, Bb, Pa, Pb, Xatol, Xbtol, Xptol);
    *MWD = calc_MWD(Mtemp, *Mn, Aa, Ab, Ba, Bb, Pa, Pb, Xatol, Xbtol, Xptol);
	*conv = 1 - ((Xa + Xb)*1.0) / (2 * VNa * conc);
    
    return;
    }

static PyObject *py_KMC(PyObject *self, PyObject *args) { // Wrap the C-based KMC function into a python-callable function
	PyObject *Py_params;	
	if (!PyArg_ParseTuple(args, "O", &Py_params)) {
		return NULL;
	}
	double params[16];
	for (int i = 0;i < 16;i++) {
		PyObject *element = PyList_GetItem(Py_params, i);
		params[i] = PyFloat_AsDouble(element);
	}

	double Mn, MWD, conv;
	KMC(params, &Mn, &MWD, &conv);
	return Py_BuildValue("ddd", Mn, MWD, conv);
}

static PyMethodDef KMCMethods[] = {

   {"KMC", py_KMC, METH_VARARGS, "simulate Mn, MWD and conversion"},
   {NULL, NULL, 0, NULL}
};

static struct PyModuleDef kmcmodule = {
	PyModuleDef_HEAD_INIT,
	"KMC",
	"KMC simulation of FRP",
	-1,
	KMCMethods
};

PyMODINIT_FUNC PyInit_KMC(void) {
	return PyModule_Create(&kmcmodule);
}

