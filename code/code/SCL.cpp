//SCL-翻转
//加载头文件
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mex.h>
//SCL译码，在此基础上加入翻转模块，即根据第一次SCL译码记录Ei值从小到大的下标，然后依次进行翻转直到通过CRC校验或者达到设置的最大翻转次数
//声明必要的函数
double f(double a,double b);//f函数
double g(double a,double b,int c);//g函数
void sort1(double *a,int *b,int len);//从小到大
void sort2(double *a,int *b,int len);//从大到小
double E1(double *arr1,double *arr2,int L);
double E2(double *arr1,double *arr2,int L);//返回值    将E值存放到一个数组中，用sort2进行排序 L表示List
int sign(double a);//delta函数 大于0判为0 小于0判为1
double **SCL(double *llr,int *CBR,int N,int n,int L);//SCL译码
/*llr表示接收似然比  CBR表示冻结比特与信息比特位置 N码长 n表示log2N L表示List 
index表示翻转的次数，比如index==0无需翻转，index==2，则需要根据Ei值从小到大进行翻转2次，即在相应的位置保留PM值相对来说较大的L条路径
对于冻结位，无需翻转，肯定正确,其中在进行第一次SCL译码时，Ei就已知了，并且我们不需要记住Ei值，只需要Ei从小到大记录相应的位置即可*/
//double *SCL(double *llr,int *CBR,int N,int n,int L,int index){
//对于C语言数组而言,下标从0开始 奇数下标执行g+mod运算 偶数下标执行g+f运算
//mod与f执行次数，与下标有关系，具体规律见下面公式，mod执行次数用数组A表示 f执行次数用数组B表示
double **SCL(double *llr,int *CBR,int N,int n,int L){
int A[N];//mod次数
int B[N];//f次数
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
double P[N][n+1][2*L];//存放中间运算似然比的值
double Q[N][n+1][2*L];//存放P的转换值，即达到2L条分枝根据PM值的大小删L条分枝需要借助的数组
int C[N][n+1][2*L];//存放bit值 
int D[N][n+1][2*L];//存放C数组转换值，原理与Q数组一样，也是进行删L条分枝时需要借助的数组
double arr1[L];
double arr2[L];
//定义一个标志flag 它的作用是记录一下什么时候达到2L条分枝。起始值flag=1 达到2L条路径后flag=2L，并且开始记录Ei的值 
//Ei的下标需要用一个数组E_arr记录
int flag=1;//是否达到2L条路径的标志 每次遇到信息位乘2 直到2L时 即flag=1
double Arr[N];//存放Ei值
double Brr[N];//存放KS值
int A_arr[N];//长度现在规定为20，即假设翻转的次数最大设定为20
int B_arr[N];
int cc=0;//的下标
//接收来的似然比llr赋值给P数组的最后1列,2L条路径都赋上 从右往左
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
//定义delta数组表示实际值，PM表示惩罚值
double del[2*L];
double PM[2*L];
//假设一开始PM值全部为0
for(int i=0;i<2*L;i++){
    PM[i]=0;
}
del[0]=P[0][0][0];
//第一位是冻结位
C[0][0][0]=0;
if(C[0][0][0]==sign(del[0])){
    PM[0]=0;
}else{
    PM[0]=fabs(del[0]);
}
//这个需要暂存值temp来存储PM值的值，因为需要根据PM值的大小进行删枝操作
double temp[2*L]; 
for(int i=0;i<flag;i++){
    temp[i]=PM[i];
}
for(int k=1;k<N;k++){
    /**********奇数操作************/
    if(k%2==1){//奇数执行g+mod
    //   printf("%d ",k);
    if(!CBR[k]){//冻结比特
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
    else{//信息位，需要展开
    flag*=2;
    if(flag<2*L){
        for(int i=0;i<flag/2;i++){
            del[i]=g(P[k-1][1][i],P[k][1][i],C[k-1][0][i]);//实际值
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
    //给Temp赋值
    for(int i=0;i<flag;i++){
        temp[i]=PM[i];//下面计算会用到
    }
  //旧条*2 给新路径赋值 P数组 C数组的值需要更新 先存放
    for(int i=0;i<k;i++){
        for(int l=0;l<flag;l++){
            D[i][0][l]=C[i][0][l/2];
        }
    }
    //再赋值回去 此时C数组已经更新完了
    for(int i=0;i<k;i++){
        for(int l=0;l<flag;l++){
            C[i][0][l]=D[i][0][l];
        }
    }
    //接着更新P数组
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
    }//大于2L
    else{
    for(int i=0;i<flag/2;i++){
        del[i]=g(P[k-1][1][i],P[k][1][i],C[k-1][0][i]);//实际值
    }
    
    for(int l=0;l<flag;l++){
        C[k][0][l]=(l%2==0)?0:1;
        if(C[k][0][l]==sign(del[l/2])){
            PM[l]=temp[l/2];
        }else{
            PM[l]=temp[l/2]+fabs(del[l/2]);
        }
    }//排序完才可以计算Ei值，下面进行删枝操作
    
    int Bra[2*L];//下标
    int Bra_1[2*L];//排序完后的下标
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
        temp[i]=PM[i];//赋值
    }//更新C数组
   
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
    }//更新P数组删枝
   
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
    else{//偶数 g+f运算
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
    if(!CBR[k]){//冻结
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
    else{//信息位
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
    }//给temp赋值
    for(int i=0;i<flag;i++){
        temp[i]=PM[i];
    }//更新C数组
    for(int i=0;i<k;i++){
        for(int l=0;l<flag;l++){
            D[i][0][l]=C[i][0][l/2];
        }
    }
    for(int i=0;i<k;i++){
        for(int l=0;l<flag;l++){
            C[i][0][l]=D[i][0][l];
        }
    }//完成C数组更新
    //接着更新P数组
    for(int l=0;l<flag;l++){
        for(int j=1;j<=n;j++){
            for(int i=0;i<N;i++){
                Q[i][j][l]=P[i][j][l/2];
            }
        }
    }//再赋值回去，P数组更新完毕
    for(int l=0;l<flag;l++){
        for(int j=1;j<=n;j++){
            for(int i=0;i<N;i++){
                P[i][j][l]=Q[i][j][l];
            }
        }
    }
    }
    else{//大于2L条路径 进行删枝
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
    int Bra[2*L];//下标
    int Bra_1[2*L];//排序完后的下标
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
        temp[i]=PM[i];//赋值
    }//更新C数组
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
    }//更新P数组删枝
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
    b[i]=(double *)malloc((L+2)*sizeof(double));//申请N*(L+2)空间 前L条路径放L条路径的译码值，最后1列数据放根据Ei值排序的下标
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
//实现函数
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
//从小到大
void sort1(double *a,int *b,int len){
    double temp1;
    int temp2;//下标
    for(int j=0;j<len;j++){
        for(int i=0;i<len-1-j;i++){
            if(a[i]>a[i+1]){
                temp1=a[i];
                a[i]=a[i+1];
                a[i+1]=temp1;
                //下标
                temp2=b[i];
                b[i]=b[i+1];
                b[i+1]=temp2;
            }
        }
    }
}
//从大到小
void sort2(double *a,int *b,int len){
    double temp1;
    int temp2;//下标
    for(int j=0;j<len;j++){
        for(int i=0;i<len-1-j;i++){
            if(a[i]<a[i+1]){
                temp1=a[i];
                a[i]=a[i+1];
                a[i+1]=temp1;
                //下标
                temp2=b[i];
                b[i]=b[i+1];
                b[i+1]=temp2;
            }
        }
    }
}
//接口函数
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]){
//与matlab中function函数类似 右边五个参数 llr CBR N n(log2(N)) L
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
