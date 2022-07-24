#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define RAND_SEED 20190719

#define N_DATA 200
#define N_TEST_DATA_SET 20

#define SIGMA 0.7 // 加わるノイズの標準偏差
#define P00 0.99 // 状態遷移確率の定義．P01 は 0→1 の遷移確率
#define P01 0.01
#define P10 0.03
#define P11 0.97

int x[N_DATA]; // もともとの信号． 0 か 1
int xmap[N_DATA]; // 推定値． 0 か 1
double y[N_DATA];  // 観測データ
int ham = 0; //hamming distance

int  Txmap[N_TEST_DATA_SET][N_DATA]; // テスト用：推定値． 
double Ty[N_TEST_DATA_SET][N_DATA];  // テスト用：観測データ

int S[N_DATA][2];  // S[2][b] = argmax_a { ( C[2][a] + h(a,b) } 
double C[N_DATA][2];

double lp[2][2]; // log(P00) などを代入．

double nrnd();

void generate_x (){

  int i;

  if ( drand48() < 0.5 ){
    x[0]=0;
  }
  else{
    x[0]=1;
  }
        
  for (i=1; i<N_DATA; i++){
    if ( x[i-1] == 0 ){
      if ( drand48() < P00 ){
 x[i] = 0;
      }
      else{
 x[i] = 1;
      }
    }
    else{
      if ( drand48() < P11 ){
 x[i] = 1;
      }
      else{
 x[i] = 0;
      }
    }
  }

}

void generate_y (){

  int i;
  for (i=0; i<N_DATA; i++){
    y[i]=(double)x[i] + SIGMA*nrnd();
  }

}


void compute_xmap (){

  int i;

  lp[0][0] = log(P00);
  lp[1][0] = log(P10);
  lp[1][1] = log(P11);
  lp[0][1] = log(P01);
    
  C[0][0] =  -(y[i]-0.0)*(y[i]-0.0)/(2.0*SIGMA*SIGMA);
  C[0][1] =  -(y[i]-1.0)*(y[i]-1.0)/(2.0*SIGMA*SIGMA);
    
  for (i=1; i<N_DATA-1; i++){
    //x_i = 0
    if(C[i-1][0] + lp[0][0] - (y[i]-0.0)*(y[i]-0.0)/(2.0*SIGMA*SIGMA) >=  C[i-1][1] + lp[1][0] - (y[i]-0.0)*(y[i]-0.0)/(2*SIGMA*SIGMA)){
      C[i][0] = C[i-1][0] + lp[0][0] - (y[i]-0.0)*(y[i]-0.0)/(2.0*SIGMA*SIGMA);
      S[i-1][0] = 0;
    }else{
      C[i][0] = C[i-1][1] + lp[1][0] - (y[i]-0.0)*(y[i]-0.0)/(2*SIGMA*SIGMA);
      S[i-1][0] = 1;
    }
      
    //x_i = 1
    if ( C[i-1][1] + lp[1][1] - (y[i]-1.0)*(y[i]-1.0)/(2.0*SIGMA*SIGMA) >=  C[i-1][0] + lp[0][1] - (y[i]-1.0)*(y[i]-1.0)/(2*SIGMA*SIGMA)){
      C[i][1] = C[i-1][1] + lp[1][1] - (y[i]-1.0)*(y[i]-1.0)/(2.0*SIGMA*SIGMA);
      S[i-1][1] = 1;
    }else{
      C[i][1] = C[i-1][0] + lp[0][1] - (y[i]-1.0)*(y[i]-1.0)/(2*SIGMA*SIGMA);
      S[i-1][1] = 0;
    }

  }
    

  if(C[N_DATA-1][0] >= C[N_DATA-1][1]){
    xmap[N_DATA-1] = 0;
    for (i=N_DATA-2; i>=0; i--){
      xmap[i] = S[i][xmap[i+1]];
    }
  }else{
    xmap[N_DATA-1] = 1;
    for (i=N_DATA-2; i>=0; i--){
      xmap[i] = S[i][xmap[i+1]];
    }
  }

  for(i=0; i<N_DATA; i++){
    xmap[i] = xmap[i];
  }

}


void ham_dis(){
  
  int i;

  for(i=0; i<N_DATA; i++){
    ham = ham + (x[i] - xmap[i])*(x[i] - xmap[i]);
  }

  printf("\nhamming_distance = %d\n",ham);
}



void show_resuls(){
    
  int i;
  for (i=0; i<N_DATA; i++){
    printf("%d\t%d\t%.8lf\t%d\n",i, x[i],y[i],xmap[i]+3 );
  }

}



/* 標準正規分布にしたがう擬似乱数の生成 */
double nrnd(){
  static int sw=0;
  static double r1,r2,s;
    
  if (sw==0){
    sw=1;
    do {
      r1=2.0*drand48()-1.0;
      r2=2.0*drand48()-1.0;
      s=r1*r1+r2*r2;
    } while (s>1.0 || s==0.0);
    s=sqrt(-2.0*log(s)/s);
    return(r1*s);
  }
  else {
    sw=0;
    return(r2*s);
  }
}


void test20(){

  int i,j;
  int n_passed;
  
  // テスト用データ 20例のファイルからの読み込み
  
  FILE *fp;
  fp = fopen("r20190702_20_test_cases", "r");
  if (fp == NULL){
    fprintf(stderr, "Can't open data file.\n");
    exit(-1);
  } 
  for (i=0; i<N_DATA; i++) {
    for (j=0; j<N_TEST_DATA_SET; j++) {
      fscanf(fp,"%lf",&Ty[j][i]);
      fscanf(fp,"%d",&Txmap[j][i]);
    }
  } 
  fclose(fp);
  
  for (j=0; j<N_TEST_DATA_SET; j++) {
    
    for (i=0; i<N_DATA; i++) {
      y[i] = Ty[j][i];
    }
    compute_xmap ();
    n_passed = 0;
    for (i=0; i<N_DATA; i++) {
      if ( xmap[i] + 3 == Txmap[j][i] ){
        n_passed++;
      }
    }
    if ( n_passed == N_DATA){
      printf("%d: passed\n", j);
    }
    else{
      printf("%d: not passed\n", j);
    }
  } 
  
}

    
void demo(){

  /* 問題を作る（200 個のデータ生成） */
  generate_x ();
  generate_y ();
  
  /* 復元する */
  compute_xmap (); 
  
  /* 結果を表示する */
  show_resuls();

  /*output ham*/
  ham_dis();

}


int main ( int argc , char * argv []){

  int i;
  int seed = RAND_SEED;
  /*
    for (i=1; i<argc; i++) {
    switch (*(argv[i]+1)) {
    case 'r':
    seed = atoi(argv[++i]);
    break;
    default:
    fprintf(stderr, "Usage : %s\n",argv[0]);
    fprintf(stderr, "\t-r : random seed(%d)\n",seed);
    exit(0);
    break;
    }
    }*/
  srand48(seed); /* 擬似乱数の種を設定 */

  //test20();
  demo();

  return 0;

}