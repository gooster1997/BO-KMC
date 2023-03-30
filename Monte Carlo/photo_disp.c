
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <unistd.h>
#include <time.h>
#include <Python.h>

void calc_c(double k[], double c[], int *V, double *conc){

    c[0]=k[0];
    c[1]=k[1]*(*conc)/(*V);
    c[2]=k[2]*(*conc)/(*V);
    c[3]=k[3]*(*conc)/(*V);
    c[4]=k[4]*(*conc)/(*V);
    c[5]=k[5]*(*conc)/(*V);
    c[6]=k[6]*(*conc)/(*V);
    c[7]=k[7]*(*conc)/(*V);
    c[8]=2.0*k[8]*(*conc)/(*V);
    c[9]=2.0*k[9]*(*conc)/(*V);
    c[10]=k[10]*(*conc)/(*V);
    c[11]=k[11]*(*conc)/(*V);
    c[12]=k[12]*(*conc)/(*V);
    c[13]=k[13]*(*conc)/(*V);

    return;
    }

void calc_a(int Xi1, int Xi2, int Xa, int Xb, long Xatol, long Xbtol, double c[], double a[], double conc, int V){
    
	double i;
	i = -(10.0)*Xi1*(conc) / V;
	a[0]=10000*c[0]*Xi1* pow(exp(1), i);
    a[1]=Xi2*c[1]*Xa;
    a[2]=Xi2*c[2]*Xb;
    a[3]=Xa*Xatol*c[3];
    a[4]=Xb*Xatol*c[4];
    a[5]=Xa*Xbtol*c[5];
    a[6]=Xb*Xbtol*c[6];
    a[7]=Xatol*Xbtol*c[7];
    a[8]=Xatol*(Xatol-1)*c[8]/2;
    a[9]=Xbtol*(Xbtol-1)*c[9]/2;
    a[10]=Xa*Xatol*c[10];
    a[11]=Xb*Xatol*c[11];
    a[12]=Xa*Xbtol*c[12];
    a[13]=Xb*Xbtol*c[13];
    a[14]=0;
    
    for (int j=0; j<=13; j++) {a[14]=a[14]+a[j];}
    
    return;
    }

int gene_mu(double a[]){
    double r;
    r=rand()/(RAND_MAX+1.0);
    double sum=0;
    int i=0;
    double temp;
    temp=r*a[14];
    while(sum<temp){
        sum=sum+a[i];
        i++;
        }
    
    return i;
    }

int chain_judge(long Xtol){
    double r=0;
    r=rand()/(RAND_MAX+1.0);
    int num;
    num=Xtol*r;

    return num;
    }

void homopro(long X[], long Xtol){
    int num;
    num=chain_judge(Xtol);
    X[num]++;

    return;
    }

int chain_delete(long X[], int num, long Xtol){
    int n;
    n=X[num];
    for (int i=num;i<=Xtol-2;i++) {X[i]=X[i+1];}
    return n;
    }

void crosspro(long Aa[], long Ab[], long Ba[], long Bb[], long Xatol, long Xbtol){
    int num, a, b;
    num=chain_judge(Xatol);
    a=chain_delete(Aa, num, Xatol);
    b=chain_delete(Ab, num, Xatol);
    Ba[Xbtol]=a;
    Bb[Xbtol]=b+1;

    return;
    }

void homoter(long a[], long b[], long Pa[], long Pb[], long Xtol, int Xptol) {
	int num1, num2, a1, a2, b1, b2;
	num1 = chain_judge(Xtol);
	a1 = chain_delete(a, num1, Xtol);
	b1 = chain_delete(b, num1, Xtol);
	num2 = chain_judge(Xtol - 1);
	a2 = chain_delete(a, num2, Xtol - 1);
	b2 = chain_delete(b, num2, Xtol - 1);
	Pa[Xptol] = a1;
	Pb[Xptol] = b1;
	Pa[Xptol + 1] = a2;
	Pb[Xptol + 1] = b2;

	return;
}

void crosster(long Aa[], long Ab[], long Ba[], long Bb[], long Pa[], long Pb[], long Xatol, long Xbtol, int Xptol) {
	int num1, num2, a1, a2, b1, b2;
	num1 = chain_judge(Xatol);
	num2 = chain_judge(Xbtol);
	a1 = chain_delete(Aa, num1, Xatol);
	b1 = chain_delete(Ab, num1, Xatol);
	a2 = chain_delete(Ba, num2, Xbtol);
	b2 = chain_delete(Bb, num2, Xbtol);
	Pa[Xptol] = a1;
	Pb[Xptol] = b1;
	Pa[Xptol + 1] = a2;
	Pb[Xptol + 1] = b2;

	return;
}

void trn(long A[], long B[], long Pa[], long Pb[], long Xtol, int Xptol){
    int num, a, b;
    num=chain_judge(Xtol);
    a=chain_delete(A, num, Xtol);
    b=chain_delete(B, num, Xtol);
    Pa[Xptol]=a;
    Pb[Xptol]=b;

    return;
    }

void reaction_execute(int mu, int *Xi1, int *Xi2, long *Xa, long *Xb, long *Xatol, long *Xbtol, int *Xptol, long Aa[], long Ab[], long Ba[], long Bb[], long Pa[], long Pb[]){

    switch(mu){
        case 1:
            (*Xi1)--;
            (*Xi2)+=2;
            break;
        case 2:
            (*Xi2)--;
            (*Xa)--;
            Aa[*Xatol]=1;
            Ab[*Xatol]=0;
            (*Xatol)++;
            break;
        case 3:
            (*Xi2)--;
            (*Xb)--;
            Bb[*Xbtol]=1;
            Ba[*Xbtol]=0;
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
			(*Xptol) += 2;
			break;
		case 9:
			homoter(Aa, Ab, Pa, Pb, *Xatol, *Xptol);
			(*Xatol) -= 2;
			(*Xptol) += 2;
			break;
		case 10:
			homoter(Ba, Bb, Pa, Pb, *Xbtol, *Xptol);
			(*Xbtol) -= 2;
			(*Xptol) += 2;
			break;
        case 11:
            trn(Aa, Ab, Pa, Pb, *Xatol, *Xptol);
            (*Xa)--;
            (*Xi2)++;
            break;
        case 12:
            trn(Aa, Ab, Pa, Pb, *Xatol, *Xptol);
            (*Xb)--;
            (*Xi2)++;
            break;
        case 13:
            trn(Ba, Bb, Pa, Pb, *Xbtol, *Xptol);
            (*Xa)--;
            (*Xi2)++;
            break;
        case 14:
            trn(Ba, Bb, Pa, Pb, *Xbtol, *Xptol);
            (*Xb)--;
            (*Xi2)++;
            break;
        }
    }

double calc_Mn(long *Mtemp, long Aa[], long Ab[], long Ba[], long Bb[], long Pa[], long Pb[], long Xatol, long Xbtol, int Xptol){
    double Mn;

    for (int i=0;i<Xatol;i++){*Mtemp=*Mtemp+150*Aa[i]+200*Ab[i];}
    for (int i=0;i<Xbtol;i++){*Mtemp=*Mtemp+150*Ba[i]+200*Bb[i];}
    for (int i=0;i<Xptol;i++){*Mtemp=*Mtemp+150*Pa[i]+200*Pb[i];}

    Mn=1.0*(*Mtemp)/(Xatol+Xbtol+Xptol);
    return Mn;
    }

double calc_MWD(long Mtemp, double Mn, long Aa[], long Ab[], long Ba[], long Bb[], long Pa[], long Pb[], long Xatol, long Xbtol, int Xptol){
    long long M2;
    double MWD, Mw;
    M2=0;

    for (int i=0;i<Xatol;i++){M2=M2+(150*Aa[i]+200*Ab[i])*(150*Aa[i]+200*Ab[i]);}
    for (int i=0;i<Xbtol;i++){M2=M2+(150*Ba[i]+200*Bb[i])*(150*Ba[i]+200*Bb[i]);}
    for (int i=0;i<Xptol;i++){M2=M2+(150*Pa[i]+200*Pb[i])*(150*Pa[i]+200*Pb[i]);}

    Mw=1.0*M2/Mtemp;
    MWD=1.0*Mw/Mn;
    return MWD;
    }

static PyObject *Monte_result(PyObject *self, PyObject *args){

    PyObject *klist;
    
    if (!PyArg_ParseTuple(args, "O", &klist)){
        return NULL;
    }
    
    double k[16];
    for (int i=0;i<16;i++){
        PyObject *temp=PyList_GetItem(klist, i);
        k[i]=PyFloat_AsDouble(temp);
    }
    
    double c[14]={0};
    double a[15]={0};
    long Aa[10000000]={0};
    long Ab[10000000]={0};
    long Ba[10000000]={0};
    long Bb[10000000]={0};
    long Pa[10000000]={0};
    long Pb[10000000]={0};
    
    int Xi1, Xi2, V, Xptol, Xi1thre;
    double conc, ratio;
	ratio = k[14];
	conc = k[15];

    Xi1=5000000*ratio;
    Xi1thre=2500000*ratio;
    Xi2=Xptol=0;
    V=5000000;
    long Xatol, Xbtol, Xa, Xb;
    Xatol=Xbtol=0;
    Xa=Xb=5000000;

    int mu;
    long cycle=1;
    long Mtemp;
    double Mn, MWD, conv;

    time_t start;
    start=time(NULL);

    srand(2021);

    calc_c(k, c, &V, &conc);

    while((Xa+Xb)>=0 && Xi1>Xi1thre){
        calc_a(Xi1, Xi2, Xa, Xb, Xatol, Xbtol,c,a,conc,V);
        mu=gene_mu(a);
        reaction_execute(mu, &Xi1, &Xi2, &Xa, &Xb, &Xatol, &Xbtol, &Xptol, Aa, Ab, Ba, Bb, Pa, Pb);
        
        if((time(NULL)-start)>900) {break;}

        cycle++;
        }
    Mtemp=0;
    Mn=calc_Mn(&Mtemp, Aa, Ab, Ba, Bb, Pa, Pb, Xatol, Xbtol, Xptol);
    MWD=calc_MWD(Mtemp, Mn, Aa, Ab, Ba, Bb, Pa, Pb, Xatol, Xbtol, Xptol);
	conv = 1 - ((Xa + Xb)*1.0) / 10000000;
    
    return Py_BuildValue("ddd", Mn, MWD, conv);
    }
    
static PyMethodDef Monte_funcs[]={

    {"photo_disp", Monte_photo_disp, METH_VARARGS, "Mn and MWD"},
    {NULL, NULL, 0, NULL}
};    

static struct PyModuleDef resultmodule={
    PyModuleDef_HEAD_INIT,
    "Monte",
    "Python interface for the Monte under C",
    -1,
    Monte_funcs
};

PyMODINIT_FUNC PyInit_Monte(void){
    return PyModule_Create(&resultmodule);
}

