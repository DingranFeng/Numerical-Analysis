#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define n 11
int i,j,s,k,ll,p;
double A[n][n],q[n][n],r[n][n],rq[n][n],I[n][n],CC[2][n],M[n][n];
double P[n],W[n],u[n],Q[n],RR[n],v[n];
double r1,r2,r3,r4,r5,s1,t,x,sum;
int n1,n2;
void Amat() {
	for(i=1; i<11; i++) {
		for(j=1; j<11; j++) {
			if(j!=i)
				A[i][j]=sin(0.5*i+0.2*j);
			else
				A[i][j]=1.5*cos(i+1.2*j);
		}
	}
}
void UpTri() {
	for(s=1; s<n-2; s++) {
		for(r4=0.0,i=s+2; i<n; i++)
			r4+=A[i][s]*A[i][s];//printf("r4=%1.12e\n",r4);
		if(r4==0)
			continue;
		else {
			r4+=A[s+1][s]*A[s+1][s];
			r1=sqrt(r4);//printf("r1=%1.12e\n",r1);
			if(A[s+1][s]>0) r2=-r1;
			else r2=r1;
			r3=r2*r2-r2*A[s+1][s];
			for(i=1; i<=s; i++)
				u[i]=0.0;
			u[s+1]=A[s+1][s]-r2;
			for(i=s+2; i<n; i++)
				u[i]=A[i][s];
			for(j=1; j<n; j++) { //P
				for(P[j]=0.0,i=1; i<n; i++)
					P[j]+=A[i][j]*u[i]/r3;
			}
			for(r5=0.0,i=1; i<n; i++) { //r5
				r5+=P[i]*u[i]/r3;
			}
			for(i=1; i<n; i++) { //Q
				for(Q[i]=0.0,j=1; j<n; j++)
					Q[i]+=A[i][j]*u[j]/r3;
			}
			for(i=1; i<n; i++) {
				W[i]=Q[i]-r5*u[i];
			}
			for(i=1; i<n; i++) {
				for(j=1; j<n; j++)
					A[i][j]-=(W[i]*u[j]+u[i]*P[j]);
			}
		}
	}
	for(i=1; i<n; i++) {
		for(j=1; j<n; j++)
			if(fabs(A[i][j])<0.000000000001)
				A[i][j]=0.0;
	}
	printf("\n进行拟上三角化后:\n");
	for(i=1; i<n; i++) {
		for(j=1; j<n; j++)
			printf("%1.12e,",A[i][j]);
		printf("\n");
	}
}
void QRgen() {
	double u[n],w[n],F[n];
	for(s=1; s<n; s++)
		for(p=1; p<n; p++)
			r[s][p]=A[s][p];
	for(s=1; s<n-1; s++) {
		for(r1=0.0,i=s; i<n; i++)
			r1+=r[i][s]*r[i][s];
		r1=sqrt(r1);
		if (r1==0)
			continue;
		else {
			if(A[s][s]>0) r2=-r1;
			else r2=r1;
			r3=r2*r2-r2*r[s][s];
			for(i=1; i<s; i++)
				u[i]=0;
			u[s]=r[s][s]-r2;
			for(i=s+1; i<n; i++)
				u[i]=r[i][s];
			for(i=1; i<n; i++) {
				for (F[i]=0.0,j=1; j<n; j++)
					F[i]+=r[j][i]*u[j]/r3;
			}
			for(i=1; i<n; i++) {
				for(j=1; j<n; j++)
					r[i][j]=r[i][j]-u[i]*F[j];
			}
			for(i=1; i<n; i++) {
				for(w[i]=0.0,j=1; j<n; j++)
					w[i]+=q[i][j]*u[j];
			}
			for(i=1; i<n; i++) {
				for(j=1; j<n; j++)
					q[i][j]-=w[i]*u[j]/r3;
			}
		}
	}
	for(i=1; i<n; i++) {
		for(j=1; j<n; j++) {
			for(rq[i][j]=0.0,s=1; s<n; s++)
				rq[i][j]+=r[i][s]*q[s][j];
		}
	}
	printf("\nQ矩阵：\n");
	for(i=1; i<n; i++) {
		for(j=1; j<n; j++)
			if(fabs(q[i][j])<0.000000000001)
				q[i][j]=0.0;
	}
	for(i=1; i<n; i++) {
		for(j=1; j<n; j++)
			printf("%1.12e,",q[i][j]);
		printf("\n");
	}
	printf("\nR矩阵：\n");
	for(i=1; i<n; i++) {
		for(j=1; j<n; j++)
			if(fabs(r[i][j])<0.000000000001)
				r[i][j]=0.0;
	}
	for(i=1; i<n; i++) {
		for(j=1; j<n; j++)
			printf("%1.12e,",r[i][j]);
		printf("\n");
	}
	printf("\nRQ矩阵：\n");
	for(i=1; i<n; i++) {
		for(j=1; j<n; j++)
			if(fabs(rq[i][j])<0.000000000001)
				rq[i][j]=0.0;
	}
	for(i=1; i<n; i++) {
		for(j=1; j<n; j++)
			printf("%1.12e,",rq[i][j]);
		printf("\n");
	}
}
void QR() {
	int K=1,m=10;
	n1=0,n2=0;
loop3:
	if(fabs(A[m][m-1])<=0.000000000001) {
		n1++;
		RR[n1]=A[m][m];
		m--;
		goto loop4;
	} else goto loop5;
loop4:
	if(m==1) {
		n1++;
		RR[n1]=A[1][1];
		goto loop9;
	} else if(m==2) {
		s1=A[m-1][m-1]+A[m][m];
		t=A[m-1][m-1]*A[m][m]-A[m][m-1]*A[m-1][m];
		x=s1*s1-4*t;
		if(x>=0) {
			n1++;
			RR[n1]=(s1+sqrt(x))/2;
			n1++;
			RR[n1]=(s1-sqrt(x))/2;
		} else {
			n2++;
			CC[0][n2]=s1/2;
			CC[1][n2]=sqrt(-x)/2;
			n2++;
			CC[0][n2]=s1/2;
			CC[1][n2]=-sqrt(-x)/2;
		}
		goto loop9;
	} else goto loop3;
loop5: {
		if(fabs(A[m-1][m-2])<=0.000000000001) {
			s1=A[m-1][m-1]+A[m][m];
			t=A[m-1][m-1]*A[m][m]-A[m][m-1]*A[m-1][m];
			x=s1*s1-4*t;
			if(x>=0) {
				n1++;
				RR[n1]=(s1+sqrt(x))/2;
				n1++;
				RR[n1]=(s1-sqrt(x))/2;
			} else {
				n2++;
				CC[0][n2]=s1/2;
				CC[1][n2]=sqrt(-x)/2;
				n2++;
				CC[0][n2]=s1/2;
				CC[1][n2]=-sqrt(-x)/2;
			}
			m--;
			m--;
			goto loop4;
		} else goto loop6;
	}
loop6: {
		if(K==2500)
			printf("未得到所有特征值！\n");
		else goto loop7;
	}
loop7: {
		s1=A[m-1][m-1]+A[m][m];
		t=A[m-1][m-1]*A[m][m]-A[m][m-1]*A[m-1][m];

		for(i=1; i<=m; i++)
			for(j=1; j<=m; j++) {
				for(sum=0.0,s=1; s<=m; s++)
					sum+=A[i][s]*A[s][j];
				M[i][j]=sum-s1*A[i][j]+t*I[i][j];
			}
		for(s=1; s<=m; s++) {
			for(r4=0.0,i=s+1; i<=m; i++)
				r4+=M[i][s]*M[i][s];
			if(r4==0)
				continue;
			else {
				r4+=M[s][s]*M[s][s];
				r1=sqrt(r4);
				if(M[s][s]>0)
					r2=-r1;
				else r2=r1;
				r3=r2*r2-r2*M[s][s];
				for(i=1; i<s; i++)
					u[i]=0.0;
				u[s]=M[s][s]-r2;
				for(i=s+1; i<=m; i++)
					u[i]=M[i][s];
				for(j=1; j<=m; j++)
					for(v[j]=0.0,i=1; i<=m; i++)
						v[j]+=M[i][j]*u[i]/r3;
				for(i=1; i<=m; i++)
					for(j=1; j<=m; j++)
						M[i][j]-=u[i]*v[j];
				for(j=1; j<=m; j++)
					for(P[j]=0.0,i=1; i<=m; i++)
						P[j]+=A[i][j]*u[i]/r3;
				for(i=1; i<=m; i++)
					for(Q[i]=0.0,j=1; j<=m; j++)
						Q[i]+=A[i][j]*u[j]/r3;
				for(r5=0.0,i=1; i<=m; i++)
					r5+=P[i]*u[i]/r3;
				for(i=1; i<=m; i++)
					W[i]=Q[i]-r5*u[i];
				for(i=1; i<n; i++)
					for(j=1; j<n; j++)
						A[i][j]-=(W[i]*u[j]+u[i]*P[j]);
			}
		}
		goto loop8;
	}
loop8: {
		k++;
		goto loop3;
	}
loop9: {
		;
	}
	printf("\nA的特征值：\n");
	for(i=1; i<=n1; i++)
		printf("%1.12e\n",RR[i]);
	for(i=1; i<=n2; i++) {
		printf("%1.12e+0*i",CC[0][i]);
		if(CC[1][i]>=0)
			printf("+*i%1.12e\n",CC[1][i]);
		else printf("*i%1.12e\n",CC[1][i]);
	}
}
void Gaussian() {
	double ch,m[n],x[n][n];
	for(p=1; p<=n1; p++) {
		Amat();
		for(i=1; i<n; i++)
			for(j=1; j<n; j++)
				A[i][j]-=RR[p]*I[i][j];
		for(k=1; k<n; k++) {
			for(ll=k,i=k; i<n; i++)
				if(A[ll][k]<A[i][k])
					ll=i;
			for(j=k; j<n; j++) {
				ch=A[ll][j];
				A[ll][j]=A[k][j];
				A[k][j]=ch;
			}
			for(i=k+1; i<n; i++) {
				m[i]=A[i][k]/A[k][k];
				for(j=k+1; j<n; j++)
					A[i][j]-=m[i]*A[k][j];
			}
		}
		x[p][n-1]=1.0;
		for(k=n-2; k>0; k--) {
			for(r4=0.0,j=k+1; j<=n; j++)
				r4+=A[k][j]*x[p][j];
		}
	}
	for(p=1; p<=n1; p++) {
		printf("\n特征值lambda = %1.12e对应的特征向量:\n",RR[p]);
		for(i=1; i<n; i++)
			printf("%1.12e\n",x[p][i]);
	}
////////////////////////////////////////////////////////////////////////
}
main() {
	for(i=1; i<n; i++) {
		for(j=1; j<n; j++) {
			if (i==j) {
				I[i][j]=1;
				q[i][j]=1;
			} else {
				I[i][j]=0;
				q[i][j]=0;
			}
		}
	}
	Amat();
	UpTri();
	QRgen();
	QR();
	Gaussian();
}
