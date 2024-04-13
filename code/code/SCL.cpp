//SCL-��ת
//����ͷ�ļ�
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mex.h>
//SCL���룬�ڴ˻����ϼ��뷭תģ�飬�����ݵ�һ��SCL�����¼Eiֵ��С������±꣬Ȼ�����ν��з�תֱ��ͨ��CRCУ����ߴﵽ���õ����ת����
//������Ҫ�ĺ���
double f(double a,double b);//f����
double g(double a,double b,int c);//g����
void sort1(double *a,int *b,int len);//��С����
void sort2(double *a,int *b,int len);//�Ӵ�С
double E1(double *arr1,double *arr2,int L);
double E2(double *arr1,double *arr2,int L);//����ֵ    ��Eֵ��ŵ�һ�������У���sort2�������� L��ʾList
int sign(double a);//delta���� ����0��Ϊ0 С��0��Ϊ1
double **SCL(double *llr,int *CBR,int N,int n,int L);//SCL����
/*llr��ʾ������Ȼ��  CBR��ʾ�����������Ϣ����λ�� N�볤 n��ʾlog2N L��ʾList 
index��ʾ��ת�Ĵ���������index==0���跭ת��index==2������Ҫ����Eiֵ��С������з�ת2�Σ�������Ӧ��λ�ñ���PMֵ�����˵�ϴ��L��·��
���ڶ���λ�����跭ת���϶���ȷ,�����ڽ��е�һ��SCL����ʱ��Ei����֪�ˣ��������ǲ���Ҫ��סEiֵ��ֻ��ҪEi��С�����¼��Ӧ��λ�ü���*/
//double *SCL(double *llr,int *CBR,int N,int n,int L,int index){
//����C�����������,�±��0��ʼ �����±�ִ��g+mod���� ż���±�ִ��g+f����
//mod��fִ�д��������±��й�ϵ��������ɼ����湫ʽ��modִ�д���������A��ʾ fִ�д���������B��ʾ
double **SCL(double *llr,int *CBR,int N,int n,int L){
int A[N];//mod����
int B[N];//f����
for(int i=0;i<N;i++){//mod
    int j=i+1;
    int count1=0;
    while (j%2==0){
        j/=2;
        count1++;
    }
    A[i]=count1;  
}//f
B[0]=0;
for(int j=1;j<N;j++){
    int count2=0;
    int i=j;
    while (i%2==0){
        i/=2;
        count2++;
    }
    B[j]=count2;
}
double P[N][n+1][2*L];//����м�������Ȼ�ȵ�ֵ
double Q[N][n+1][2*L];//���P��ת��ֵ�����ﵽ2L����֦����PMֵ�Ĵ�СɾL����֦��Ҫ����������
int C[N][n+1][2*L];//���bitֵ 
int D[N][n+1][2*L];//���C����ת��ֵ��ԭ����Q����һ����Ҳ�ǽ���ɾL����֦ʱ��Ҫ����������
double arr1[L];
double arr2[L];
//����һ����־flag ���������Ǽ�¼һ��ʲôʱ��ﵽ2L����֦����ʼֵflag=1 �ﵽ2L��·����flag=2L�����ҿ�ʼ��¼Ei��ֵ 
//Ei���±���Ҫ��һ������E_arr��¼
int flag=1;//�Ƿ�ﵽ2L��·���ı�־ ÿ��������Ϣλ��2 ֱ��2Lʱ ��flag=1
double Arr[N];//���Eiֵ
double Brr[N];//���KSֵ
int A_arr[N];//�������ڹ涨Ϊ20�������跭ת�Ĵ�������趨Ϊ20
int B_arr[N];
int cc=0;//���±�
//����������Ȼ��llr��ֵ��P��������1��,2L��·�������� ��������
for(int j=0;j<2*L;j++){
    for(int i=0;i<N;i++){
        P[i][n][j]=llr[i];
    }
}
for(int l=0;l<2*L;l++){
int h=N;
for(int j=n-1;j>=0;j--){
    for(int i=0;i<h/2;i++){
        P[i][j][l]=f(P[i][j+1][l],P[i+h/2][j+1][l]);
    }
    h/=2;
}    
}
//����delta�����ʾʵ��ֵ��PM��ʾ�ͷ�ֵ
double del[2*L];
double PM[2*L];
//����һ��ʼPMֵȫ��Ϊ0
for(int i=0;i<2*L;i++){
    PM[i]=0;
}
del[0]=P[0][0][0];
//��һλ�Ƕ���λ
C[0][0][0]=0;
if(C[0][0][0]==sign(del[0])){
    PM[0]=0;
}else{
    PM[0]=fabs(del[0]);
}
//�����Ҫ�ݴ�ֵtemp���洢PMֵ��ֵ����Ϊ��Ҫ����PMֵ�Ĵ�С����ɾ֦����
double temp[2*L]; 
for(int i=0;i<flag;i++){
    temp[i]=PM[i];
}
for(int k=1;k<N;k++){
    /**********��������************/
    if(k%2==1){//����ִ��g+mod
    //   printf("%d ",k);
    if(!CBR[k]){//�������
    for(int i=0;i<flag;i++){
        del[i]=g(P[k-1][1][i],P[k][1][i],C[k-1][0][i]);
        C[k][0][i]=0;
        if(C[k][0][i]==sign(del[i])){
            PM[i]=temp[i];
        }else{
            PM[i]=temp[i]+fabs(del[i]);
        }
    }
    for(int i=0;i<flag;i++){
        temp[i]=PM[i];
    }
    }
    else{//��Ϣλ����Ҫչ��
    flag*=2;
    if(flag<2*L){
        for(int i=0;i<flag/2;i++){
            del[i]=g(P[k-1][1][i],P[k][1][i],C[k-1][0][i]);//ʵ��ֵ
        }
    for(int l=0;l<flag;l++){
        C[k][0][l]=(l%2==0)?0:1;
        if(C[k][0][l]==sign(del[l/2])){
            PM[l]=temp[l/2];
        }   
        else{
            PM[l]=temp[l/2]+fabs(del[l/2]);
        }
    }
    //��Temp��ֵ
    for(int i=0;i<flag;i++){
        temp[i]=PM[i];//���������õ�
    }
  //����*2 ����·����ֵ P���� C�����ֵ��Ҫ���� �ȴ��
    for(int i=0;i<k;i++){
        for(int l=0;l<flag;l++){
            D[i][0][l]=C[i][0][l/2];
        }
    }
    //�ٸ�ֵ��ȥ ��ʱC�����Ѿ���������
    for(int i=0;i<k;i++){
        for(int l=0;l<flag;l++){
            C[i][0][l]=D[i][0][l];
        }
    }
    //���Ÿ���P����
    for(int l=0;l<flag;l++){
        for(int j=1;j<=n;j++){
                for(int i=0;i<N;i++){
                    Q[i][j][l]=P[i][j][l/2];
                }
            }
        }
    for(int l=0;l<flag;l++){
        for(int j=1;j<=n;j++){
            for(int i=0;i<N;i++){
                    P[i][j][l]=Q[i][j][l];
                }
            }
        }
    }//����2L
    else{
    for(int i=0;i<flag/2;i++){
        del[i]=g(P[k-1][1][i],P[k][1][i],C[k-1][0][i]);//ʵ��ֵ
    }
    
    for(int l=0;l<flag;l++){
        C[k][0][l]=(l%2==0)?0:1;
        if(C[k][0][l]==sign(del[l/2])){
            PM[l]=temp[l/2];
        }else{
            PM[l]=temp[l/2]+fabs(del[l/2]);
        }
    }//������ſ��Լ���Eiֵ���������ɾ֦����
    
    int Bra[2*L];//�±�
    int Bra_1[2*L];//���������±�
    for(int i=0;i<2*L;i++){
        Bra[i]=i;
    }
    sort1(PM,Bra,2*L);
    for(int i=0;i<flag;i++){
        Bra_1[i]=Bra[i];
    }
    
   for(int i=0;i<L;i++){
        arr1[i]=PM[i];
        arr2[i]=PM[i+L];
    }
   
    Arr[cc]=E2(arr1,arr2,L);//2
	Brr[cc]=E1(arr1,arr2,L);//1
    A_arr[cc]=k;
	B_arr[cc]=k;
    cc++;
    for(int i=0;i<L;i++){
        temp[i]=PM[i];//��ֵ
    }//����C����
   
    for(int i=0;i<k;i++){
        for(int l=0;l<L;l++){
            D[i][0][l]=C[i][0][Bra_1[l]/2];
        }
    }
   
    for(int i=0;i<k;i++){
        for(int l=0;l<L;l++){
            C[i][0][l]=D[i][0][l];
        }
    } 
    
    for(int i=0;i<L;i++){
        D[k][0][i]=C[k][0][Bra[i]];
    }
    
    for(int i=0;i<L;i++){
        C[k][0][i]=D[k][0][i];
    }//����P����ɾ֦
   
    for(int l=0;l<flag/2;l++){
        for(int j=1;j<=n;j++){
            for(int i=0;i<=N;i++){
                Q[i][j][l]=P[i][j][Bra[l]/2];
            }
        }
    }

    for(int l=0;l<flag/2;l++){
        for(int j=1;j<=n;j++){
            for(int i=0;i<N;i++){
                    P[i][j][l]=Q[i][j][l];
                }
            }
        } 
    flag/=2;
    } //
    } 
    int n1=pow(2,A[k]);
    int n2=n1;
    for(int j=0;j<A[k];j++){
        for(int e=k-n1+1;e<=k;e=e+(2*n1)/n2){
            for(int i=e;i<e+n1/n2;i++){
                for(int l=0;l<flag;l++){
                    C[i][j+1][l]=(C[i][j][l]+C[i+n1/n2][j][l])%2;
                    C[i+n1/n2][j+1][l]=C[i+n1/n2][j][l];
                }
            }
        }
        n2/=2;
    }
    
    }
    else{//ż�� g+f����
    //printf("%d ",k);
    int num3=B[k];
    int j=num3;
    int num4=pow(2,num3);
    for(int r=0;r<num4;r++){
        for(int l=0;l<flag;l++){
        P[k+r][j][l]=g(P[k-num4+r][j+1][l],P[k+r][j+1][l],C[k-num4+r][j][l]);
        }
    }
    for(int y=0;y<flag;y++){
        for(int q=j;q>0;q--){
            int num5=pow(2,q-1);
            for(int i=0;i<num5;i++){
               P[k+i][q-1][y]=f(P[k+i][q][y], P[k + num5+i][q][y]);
            }  
        }
    }
    if(!CBR[k]){//����
        for(int l=0;l<flag;l++){
            C[k][0][l]=0;
            if(C[k][0][l]==sign(P[k][0][l])){
                PM[l]=temp[l];
            }
            else{
                PM[l]=temp[l]+fabs(P[k][0][l]);
            }
        }
        for(int i=0;i<flag;i++){
            temp[i]=PM[i];
        }
    }
    else{//��Ϣλ
    flag*=2;
    //printf("%d ",flag);
    if(flag<2*L){
    for(int i=0;i<flag/2;i++){
        del[i]=P[k][0][i];
    }
    for(int l=0;l<flag;l++){
        C[k][0][l]=((l%2==0)?0:1);
        if(C[k][0][l]==sign(del[l/2])){
            PM[l]=temp[l/2];
        }
        else{
            PM[l]=temp[l/2]+fabs(del[l/2]);
        }
    }//��temp��ֵ
    for(int i=0;i<flag;i++){
        temp[i]=PM[i];
    }//����C����
    for(int i=0;i<k;i++){
        for(int l=0;l<flag;l++){
            D[i][0][l]=C[i][0][l/2];
        }
    }
    for(int i=0;i<k;i++){
        for(int l=0;l<flag;l++){
            C[i][0][l]=D[i][0][l];
        }
    }//���C�������
    //���Ÿ���P����
    for(int l=0;l<flag;l++){
        for(int j=1;j<=n;j++){
            for(int i=0;i<N;i++){
                Q[i][j][l]=P[i][j][l/2];
            }
        }
    }//�ٸ�ֵ��ȥ��P����������
    for(int l=0;l<flag;l++){
        for(int j=1;j<=n;j++){
            for(int i=0;i<N;i++){
                P[i][j][l]=Q[i][j][l];
            }
        }
    }
    }
    else{//����2L��·�� ����ɾ֦
    for(int i=0;i<flag/2;i++){
        del[i]=P[k][0][i];
    }
    for(int l=0;l<flag;l++){
        C[k][0][l]=((l%2==0)?0:1);
        if(C[k][0][l]==sign(del[l/2])){
            PM[l]=temp[l/2];
        }
        else{
            PM[l]=temp[l/2]+fabs(del[l/2]);
        }
    } 
    int Bra[2*L];//�±�
    int Bra_1[2*L];//���������±�
    for(int i=0;i<2*L;i++){
        Bra[i]=i;
    }
    sort1(PM,Bra,2*L);
    for(int i=0;i<flag;i++){
        Bra_1[i]=Bra[i];
    }
   for(int i=0;i<L;i++){
        arr1[i]=PM[i];
        arr2[i]=PM[i+L];
    }
    Arr[cc]=E2(arr1,arr2,L);//1.5
	Brr[cc]=E1(arr1,arr2,L);//1
    A_arr[cc]=k;
	B_arr[cc]=k;
    cc++;
    for(int i=0;i<L;i++){
        temp[i]=PM[i];//��ֵ
    }//����C����
    for(int i=0;i<k;i++){
        for(int l=0;l<L;l++){
            D[i][0][l]=C[i][0][Bra_1[l]/2];
        }
    }
    for(int i=0;i<k;i++){
        for(int l=0;l<L;l++){
            C[i][0][l]=D[i][0][l];
        }
    }
    for(int i=0;i<L;i++){
        D[k][0][i]=C[k][0][Bra[i]];
    }
    for(int i=0;i<L;i++){
        C[k][0][i]=D[k][0][i];
    }//����P����ɾ֦
    for(int l=0;l<flag/2;l++){
        for(int j=1;j<=n;j++){
            for(int i=0;i<=N;i++){
                Q[i][j][l]=P[i][j][Bra[l]/2];
            }
        }
    }
    for(int l=0;l<flag/2;l++){
        for(int j=1;j<=n;j++){
            for(int i=0;i<N;i++){
                    P[i][j][l]=Q[i][j][l];
                }
            }
        } 
    flag/=2;
    }
    }
}
}
int PM_arr[cc];
int BM_arr[cc];
for(int i=0;i<cc;i++){
    PM_arr[i]=i;
	BM_arr[i]=i;
}
sort1(Arr,PM_arr,cc);
sort1(Brr,BM_arr,cc);
double **b=(double **)malloc(N*sizeof(double *));
for(int i=0;i<N;i++){
    b[i]=(double *)malloc((L+2)*sizeof(double));//����N*(L+2)�ռ� ǰL��·����L��·��������ֵ�����1�����ݷŸ���Eiֵ������±�
}
for(int j=0;j<L;j++){
    for(int i=0;i<N;i++){
        b[i][j]=C[i][0][j];
    }
}
for(int i=0;i<cc;i++){
    b[i][L]=A_arr[PM_arr[i]];
	b[i][L+1]=B_arr[BM_arr[i]];
}
for(int i=cc;i<N;i++){
	b[i][L]=0;
	b[i][L+1]=0;
}
return b;
}
//ʵ�ֺ���
double f(double a,double b){
    return a*b>0?fmin(fabs(a),fabs(b)):-fmin(fabs(a),fabs(b));
}
double g(double a,double b,int c){
    return c>0?(b-a):(b+a);
}
int sign(double a){
    int b;
    return b=a>0?0:1;
}
double E1(double *arr1,double *arr2,int len){
    double sum1=0;
    double sum2=0;
    double sum=0;
    for(int i=0;i<len;i++){
        sum1+=exp(-arr1[i]);
        sum2+=exp(-arr2[i]);
    }
    sum=log(sum1/pow(sum2,1.5));
    return sum;
}
double E2(double *arr1,double *arr2,int len){
    double sum1=0;
    double sum2=0;
    double sum=0;
    for(int i=0;i<len;i++){
        sum1+=exp(-arr1[i]);
        sum2+=exp(-arr2[i]);
    }
    sum=log(sum1/pow(sum2,1));
    return sum;
}
//��С����
void sort1(double *a,int *b,int len){
    double temp1;
    int temp2;//�±�
    for(int j=0;j<len;j++){
        for(int i=0;i<len-1-j;i++){
            if(a[i]>a[i+1]){
                temp1=a[i];
                a[i]=a[i+1];
                a[i+1]=temp1;
                //�±�
                temp2=b[i];
                b[i]=b[i+1];
                b[i+1]=temp2;
            }
        }
    }
}
//�Ӵ�С
void sort2(double *a,int *b,int len){
    double temp1;
    int temp2;//�±�
    for(int j=0;j<len;j++){
        for(int i=0;i<len-1-j;i++){
            if(a[i]<a[i+1]){
                temp1=a[i];
                a[i]=a[i+1];
                a[i+1]=temp1;
                //�±�
                temp2=b[i];
                b[i]=b[i+1];
                b[i+1]=temp2;
            }
        }
    }
}
//�ӿں���
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]){
//��matlab��function�������� �ұ�������� llr CBR N n(log2(N)) L
    double *a=mxGetPr(prhs[0]);
    double *CB=mxGetPr(prhs[1]);
    int N=mxGetScalar(prhs[2]);
    int n=mxGetScalar(prhs[3]);
    int L=mxGetScalar(prhs[4]);

    double *llr=(double*)malloc(sizeof(double)*N);
    int *CBR=(int*)malloc(sizeof(int)*N);
    for(int i=0;i<N;i++){
        llr[i]=a[i];
        CBR[i]=(int)CB[i];
    }
    double **u_d;
    u_d=SCL(llr,CBR,N,n,L);
    free(llr);
    free(CBR);
    plhs[0]=mxCreateDoubleMatrix(N,L+2, mxREAL);
    double *u=mxGetPr(plhs[0]);
    for(int i=0;i<N;i++){
        for(int j=0;j<L+2;j++){
            u[j*N+i]=u_d[i][j];
        }
    }
    free(u_d);
}
