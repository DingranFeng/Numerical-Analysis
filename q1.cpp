#include <stdio.h>
#include <math.h>
#include <stdlib.h>
double mul(double x[], double y[]) {
	double e=0;
	for (int i = 2; i < 503; ++i) {
		e=e+x[i]*y[i];
	}
	return(e);
}
double dev(double lam0, double lam1) {
	double e, d, l;
	e=fabs(lam1-lam0);
	d=fabs(lam1);
	l=e/d;
	return(l);
}
void ite1(double a[], double b, double c, double *x0, double *x1) {
	for (int i = 0; i < 501; ++i) {
		x1[i+2]=c*x0[i]+b*x0[i+1]+a[i]*x0[i+2]+b*x0[i+3]+c*x0[i+4];
	}
}
void ite2(double (*L)[501], double (*U)[3], double *x0, double *x1) {
	double xx[501];
	for (int i = 0; i < 501; ++i) {
		xx[i]=x0[i];
		for (int t = (0<(i-2)?(i-2):0); t <= i-1; ++t) {
			xx[i]=xx[i]-L[i-t][t]*xx[t];
		}
	}
	x1[500]=xx[500]/U[500][0];
	for (int i = 499; i >=0; --i) {
		x1[i]=xx[i];
		int z=((i+2)<500?(i+2):500);
		for (int t = i+1; t <= ((i+2)<500?(i+2):500); ++t) {
			x1[i]=x1[i]-U[i][t-i]*x1[t];
		}
		x1[i]=x1[i]/U[i][0];
	}
}
void LU(double a[],double b, double c, double (*L)[501] , double (*U)[3]) {
	U[0][0]=a[0];
	U[0][1]=b;
	U[0][2]=c;
	L[1][0]=b/U[0][0];
	L[2][0]=c/U[0][0];
	U[1][0]=a[1]-L[1][0]*U[0][1];
	U[1][1]=b-L[1][0]*U[0][2];
	U[1][2]=c;
	L[1][1]=(b-L[2][0]*U[0][1])/U[1][0];
	L[2][1]=c/U[1][0];
	for (int k = 2; k < 499; ++k) {
		U[k][0]=a[k]-L[2][k-2]*U[k-2][2]-L[1][k-1]*U[k-1][1];
		U[k][1]=b-L[1][k-1]*U[k-1][2];
		U[k][2]=c;
		L[1][k]=(b-L[2][k-1]*U[k-1][1])/U[k][0];
		L[2][k]=c/U[k][0];
	}
	int k=499;
	U[k][0]=a[k]-L[2][k-2]*U[k-2][2]-L[1][k-1]*U[k-1][1];
	U[k][1]=b-L[1][k-1]*U[k-1][2];
	L[1][k]=(b-L[2][k-1]*U[k-1][1])/U[k][0];
	k=500;
	U[k][0]=a[k]-L[2][k-2]*U[k-2][2]-L[1][k-1]*U[k-1][1];
}
///////////////////////////////////////////////////////////////
void PowerMethod() {
	int i,k1,k2;
	double a[501]= {0},x[505]= {0.0}, y[505]= {0.0}, b=0.16, c=-0.064, lam=0, lam0=0,lam1=0,e1=0, d=0,*x1=y, *x0=x, *t;
	for (int i = 1; i < 502; ++i) {
		a[i-1]=(1.64-i*0.024)*sin(0.2*i)-0.64*exp(0.1/i);
		x[i+1]=1;
	}
	for (i = 0; i <10000; ++i) {
		d=sqrt(mul(x0,x0));
		for(int j=2; j<503; j++) {
			x0[j]=x0[j]/d;
		}
		ite1(a,b,c,x0,x1);
		lam1=mul(x0,x1);
		e1=dev(lam0,lam1);
		if (e1<1e-12) {
			lam=lam1;
			break;
		} else
			t=x0;
		x0=x1;
		x1=t;
		lam0=lam1;
	}
	k1=i-1;
	for (int i = 0; i < 501; ++i) {
		a[i]=a[i]-lam;
		x[501+1]=i;
	}
	for (i = 0; i <100000; ++i) {
		d=sqrt(mul(x0,x0));
		for(int j=2; j<503; j++) {
			x0[j]=x0[j]/d;
		}
		ite1(a,b,c,x0,x1);
		lam1=mul(x0,x1);
		e1=dev(lam0,lam1);
		if (e1<1e-12) {
			k2=i;
			lam1=lam1+lam;
			break;
		} else {
			t=x0;
		}
		x0=x1;
		x1=t;
		lam0=lam1;
	}
	lam1=lam1+lam;
	k2=i;
	printf("利用幂法求解：\n");
	printf("第一个特征值：lam1= %10.12e\n最后一个特征值：lam501= %10.12e\n\n",lam<lam1?lam:lam1,lam>lam1?-lam:-lam1);
}
///////////////////////////////////////////////////////////////
void InversePowerMethod() {
	double a[501]= {0}, aa[501], x[501]= {0.0}, y[501]= {0.0}, l[3][501]= {0},u[501][3]= {0},(*L)[501], (*U)[3], b=0.16, c=-0.064, lam0=0, lam1=0,e1=0, d=0,*x1=y, *x0=x, *t;
	for (int i = 1; i < 502; ++i) {
		a[i-1]=(1.64-i*0.024)*sin(0.2*i)-0.64*exp(0.1/i);
		x[i-1]=1;
	}
	L=l;
	U=u;
	x0=x;
	x1=y;
	LU(a,b,c,L,U);
	printf("利用反幂法求解：\n");
	for (int i = 0; i <1000; ++i) {
		d=sqrt(mul(x0,x0));
		for(int j=0; j<501; j++) {
			x0[j]=x0[j]/d;
		}
		ite2(L,U,x0,x1);
		lam1=mul(x0,x1);
		e1=dev(lam0,lam1);
		if (e1<1e-12) {
			printf("lams = %.12e\n",1/lam1);
			break;
		} else
			t=x0;
		x0=x1;
		x1=t;
		lam0=lam1;
	}
	double detA=1, condA, lams=1/lam1;
	for (int i = 0; i < 501; ++i) {
		detA=detA*U[i][0];
	}
	printf("行列式：detA = %.12e\n",detA );
	double lam= -10.700113615016, lam501= 9.724634099619,delta, miu;
	condA=fabs(lam/lams);
	printf("条件数：condA = %.12e\n\n",condA );
	delta=(lam501-lam)/40;
	for (int k = 1; k < 40; ++k) {
		miu=lam+delta*k;
		for (int i = 0; i < 501; ++i) {
			aa[i]=a[i]-miu;
			x[i]=1;
			y[i]=2;
		}
		LU(aa,b,c,L,U);
		lam0=0;
		for (int i = 0; i <1000; ++i) {
			d=sqrt(mul(x0,x0));
			for(int j=0; j<501; j++) {
				x0[j]=x0[j]/d;
			}
			ite2(L,U,x0,x1);
			lam1=mul(x0,x1);
			e1=dev(lam0,lam1);
			if (e1<1e-12) {
				printf("lami_%d = %.12e\n",k,1/lam1+miu);
				break;
			} else
				t=x0;
			x0=x1;
			x1=t;
			lam0=lam1;
		}
	}
}
///////////////////////////////////////////////////////////////
main() {
	PowerMethod();
	InversePowerMethod();
}
