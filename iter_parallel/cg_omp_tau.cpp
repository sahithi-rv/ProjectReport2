#include<bits/stdc++.h>
#include <chrono>
#include <omp.h>
#include <TAU.h>
#define PB push_back
#define MP make_pair
#define F first
#define S second
#define RI(a) scanf("%d",&a);
#define SZ(a) (int)(a.size())
#define SET(a,b) memset(a,b,sizeof(a))
#define TR(a,t) for(__typeof(a.begin()) t=a.begin();t!=a.end();t++)
#define rep(i,l,h) for(int i=(l); i<=(h);i++)
#define repd(i,h,l) for(int i=(h);i>=(l);i--)
#define ALL(a) a.begin(),a.end()
#define DRT()  int t; cin>>t; while(t--)
#define PRSNT(a,e) (a.find(e) != a.end())
#define MINH priority_queue<int, vector<int>, greater<int> >
#define N 100009
#define MOD 1000000007
#define output1(x) cout << #x << " " << x << endl;
#define output2(x,y) cout <<x<<" "<<y <<endl;
#define inp1(x) cin >> x;
#define inp2(x,y) cin >> x >> y;
#define MAX(a,b) a>b?a:b
#define MIN(a,b) a<b?a:b
#define THREADS 4

typedef long long int  ll;

using namespace std;
using namespace std::chrono;
int n;

float EPS = 0.0001;
void get_vectorB(char* fileB, float **b){
	ifstream myfile(fileB);
	rep(i,0,n-1){
		rep(j,0,n-1){
			myfile >> b[i][j];
		}
	}
	myfile.close();

}

void partial_matmul(float **b, float **out, int tid){
	// threads = 4
	int sq = sqrt(THREADS);
	int k = n/sq;
	//cout<<"k:" << k << endl;
	//cout << "partial_matmul" << tid << endl;
	if(tid == 0){
			//cout << "partial_matmul" << tid << endl;
		for(int i = 0;i<k;i++){
			for(int j = 0;j<k;j++){
				if(i==0 && j==0){
					out[i][j] = 4*b[i][j] - b[i][j+1] - b[i+1][j];
				}else if (i==0){
					out[i][j] = 4*b[i][j] - b[i+1][j] - b[i][j-1] - b[i][j+1];
				}else if (j==0){
					out[i][j] = 4*b[i][j] - b[i+1][j] - b[i-1][j] - b[i][j+1];
				}else{
					out[i][j] = 4*b[i][j] - b[i+1][j] - b[i-1][j] - b[i][j+1] - b[i][j-1];
				}
			}
		}
	}else if(tid == 1){
			//cout << "partial_matmul" << tid << endl;
		for(int i = 0;i<k;i++){
			for(int j = 1;j<=k;j++){
				if (i==0 && j==(k)){
					out[i][j-1] = 4*b[i][j] - b[i][j-1] - b[i+1][j];
				}else if(i==0){
					out[i][j-1] = 4*b[i][j] - b[i+1][j] - b[i][j-1] - b[i][j+1];
				}else if (j==(k)){
					out[i][j-1] = 4*b[i][j] - b[i-1][j] - b[i+1][j] - b[i][j-1];
				}else{
					out[i][j-1] = 4*b[i][j] - b[i+1][j] - b[i-1][j] - b[i][j+1] - b[i][j-1];
				}
			}
		}

	}else if(tid == 2){
			//cout << "partial_matmul" << tid << endl;

		for(int i = 1;i<=k;i++){
			for(int j = 0;j<k;j++){
				if (i==(k) && j==0){
					out[i-1][j] = 4*b[i][j] - b[i][j+1] - b[i-1][j];
				}else if (j==0){
					out[i-1][j] = 4*b[i][j] - b[i+1][j] - b[i-1][j] - b[i][j+1];
				}else if(i==k){
					out[i-1][j] = 4*b[i][j] - b[i][j-1] - b[i][j+1] - b[i-1][j];
				}
				else{
					out[i-1][j] = 4*b[i][j] - b[i+1][j] - b[i-1][j] - b[i][j+1] - b[i][j-1];
				}
			}
		}
	}else{
					//cout << "partial_matmul" << tid << endl;

		for(int i = 1;i<=k;i++){
			for(int j = 1;j<=k;j++){
				if(i==k && j==k){
					out[i-1][j-1] = 4*b[i][j] - b[i][j-1] - b[i-1][j];
				}else if(i==k){
					out[i-1][j-1] = 4*b[i][j] - b[i][j-1] - b[i][j+1] - b[i-1][j];
				}else if(j==k){
					out[i-1][j-1] = 4*b[i][j] - b[i-1][j] - b[i+1][j] - b[i][j-1];		
				}else{
					out[i-1][j-1] = 4*b[i][j] - b[i+1][j] - b[i-1][j] - b[i][j+1] - b[i][j-1];					
				}
			}
		}
	}
}

void matmul(float **b, float **out){
	int i,j;
	#pragma omp parallel for schedule(static) collapse(2) private(i,j) shared(b,n, out) num_threads(THREADS) 
	for(i=0; i<=n-1;i++){
		for(j=0; j<=n-1;j++){
			if(i==0 && j==0){
				out[i][j] = 4*b[i][j] - b[i][j+1] - b[i+1][j];
			}else if (i==0 && j==(n-1)){
				out[i][j] = 4*b[i][j] - b[i][j-1] - b[i+1][j];
			}else if (i==(n-1) && j==0){
				out[i][j] = 4*b[i][j] - b[i][j+1] - b[i-1][j];
			}else if (i==(n-1) && j==(n-1)){
				out[i][j] = 4*b[i][j] - b[i][j-1] - b[i-1][j];
			}else if (i==0){
				out[i][j] = 4*b[i][j] - b[i+1][j] - b[i][j-1] - b[i][j+1];
			}else if (j==0){
				out[i][j] = 4*b[i][j] - b[i+1][j] - b[i-1][j] - b[i][j+1];
			}else if (i==(n-1)){
				out[i][j] = 4*b[i][j] - b[i][j-1] - b[i][j+1] - b[i-1][j];
			}else if (j==(n-1)){
				out[i][j] = 4*b[i][j] - b[i-1][j] - b[i+1][j] - b[i][j-1];
			}else{
				out[i][j] = 4*b[i][j] - b[i+1][j] - b[i-1][j] - b[i][j+1] - b[i][j-1];
			}
		}
	}
}

float dot_product2(float **l, float **m ){

	float out = 0.0;
	int i,j;
	#pragma omp parallel for reduction(+: out) collapse(2) shared(l,m,n) private(i,j) num_threads(THREADS) 
	for(i = 0;i<n;i++){
		//#pragma omp parallel for  shared(i,n) num_threads(THREADS)
		for(j = 0;j<n;j++){
			out += l[i][j]*m[i][j];
		}
	}
	return out;
}

float dot_product1(float **l ){

	float out = 0.0;
	int i,j;
	#pragma omp parallel for reduction(+: out) collapse(2) shared(l,n) private(i,j) num_threads(THREADS) 
	for(i = 0;i<n;i++){
		//#pragma omp parallel for  shared(i,n) num_threads(THREADS)
		for(j = 0;j<n;j++){
			out += l[i][j]*l[i][j];
		}
	}
	return out;
}

float vector_add1(float **x, float alpha, float **p){
	int i,j;
	#pragma omp parallel for collapse(2) shared(x,p,alpha,n) private(i,j) num_threads(THREADS) 
	for(i = 0;i<n;i++){
		//#pragma omp parallel for schedule(static) shared(i,n) num_threads(THREADS)
		for(j = 0;j<n;j++){
			x[i][j] = x[i][j] + alpha*p[i][j];
		}
	}
}

float vector_sub(float **r, float alpha, float **z){
	int i,j;
	#pragma omp parallel for collapse(2) shared(r,alpha,z,n) private(i,j) num_threads(THREADS) 
	for(i = 0;i<n;i++){
		//#pragma omp parallel for schedule(static) shared(i,n) num_threads(THREADS)
		for(j = 0;j<n;j++){
			r[i][j] = r[i][j] - alpha*z[i][j];
		}
	}
}

float vector_add2(float **r, float beta, float **p){
	int i,j;
	#pragma omp parallel for collapse(2) shared(r,beta,p,n) private(i,j) num_threads(THREADS) 
	for(i = 0;i<n;i++){
		//#pragma omp parallel for schedule(static) shared(i,n) num_threads(THREADS)
		for(j = 0;j<n;j++){
			p[i][j] = r[i][j] + beta*p[i][j];
		}
	}
}

int main(int argc, char *argv[]){
	if(argc<2){
		printf("Enter  arguments: input file for vector B,  Grid size\n");
		return 0;		
	}
	n = atoi(argv[2]);
	cout << n << endl;
	//float r[n][n], b[n][n], z[n][n], p[n][n], x[n][n];
	float **r, **b, **z, **p, **x;
	
	b = (float **)malloc(sizeof(float*)*n);
	rep(i,0,n-1){
		b[i] = (float *)malloc((sizeof(float))*n);
	}
	r = (float **)malloc(sizeof(float*)*n);
	rep(i,0,n-1){
		r[i] = (float *)malloc((sizeof(float))*n);
	}
	z = (float **)malloc(sizeof(float*)*n);
	rep(i,0,n-1){
		z[i] = (float *)malloc((sizeof(float))*n);
	}
	p = (float **)malloc(sizeof(float*)*n);
	rep(i,0,n-1){
		p[i] = (float *)malloc((sizeof(float))*n);
	}	
	x = (float **)malloc(sizeof(float*)*n);
	rep(i,0,n-1){
		x[i] = (float *)malloc((sizeof(float))*n);
	}
	
	get_vectorB(argv[1], b);
	float del_new = 0,  del_old;
	//init b,p, r,x;
	srand((unsigned)time(0)); 

	
	for(int i = 0;i<n;i++){
   		for(int j = 0;j<n;j++){
   			//x[i][j] = 0.1*(rand()%10)+1;
   			x[i][j] = 0;
   		}
   	}

   	high_resolution_clock::time_point t1 = high_resolution_clock::now();

   	//cout << "lleling" << endl;
   //	#pragma omp parallel for schedule(static) num_threads(THREADS)
	for(int i = 0;i<n;i++){
	//	#pragma omp parallel for schedule(static) shared(i,n) num_threads(THREADS)
   		for(int j = 0;j<n;j++){
   			r[i][j] = b[i][j];
   			p[i][j] = b[i][j];
   			
        	del_new += b[i][j]*b[i][j];
        }
    }
    #pragma omp barrier

    del_old = del_new;
	float float1 = del_new;
	float float2 = del_old;
	float float3 = sqrt(abs(float1));
	float put;
    //iterations
    float vv, alpha, beta;
    float leps = log10(EPS);
	put = log10(float3/sqrt(float2) ) ;
	int sq = sqrt(THREADS);
	int row_sz = n/sq;
    int col_sz = n/sq;
    int iter = 0;
	
	#pragma omp barrier
    
    while (put >= leps){

		int i,j;
    	#pragma omp parallel private(i,j) shared(p,n, x,r, alpha, beta, z, del_new, del_old, vv) num_threads(THREADS) 
    	{
    		TAU_PROFILE_TIMER(timer, "matmul", "", TAU_USER);
    		TAU_PROFILE_START(timer)
    		#pragma omp for schedule(static) collapse(2)  
			for(i=0; i<=n-1;i++){
				for(j=0; j<=n-1;j++){
					if(i==0 && j==0){
						z[i][j] = 4*p[i][j] - p[i][j+1] - p[i+1][j];
					}else if (i==0 && j==(n-1)){
						z[i][j] = 4*p[i][j] - p[i][j-1] - p[i+1][j];
					}else if (i==(n-1) && j==0){
						z[i][j] = 4*p[i][j] - p[i][j+1] - p[i-1][j];
					}else if (i==(n-1) && j==(n-1)){
						z[i][j] = 4*p[i][j] - p[i][j-1] - p[i-1][j];
					}else if (i==0){
						z[i][j] = 4*p[i][j] - p[i+1][j] - p[i][j-1] - p[i][j+1];
					}else if (j==0){
						z[i][j] = 4*p[i][j] - p[i+1][j] - p[i-1][j] - p[i][j+1];
					}else if (i==(n-1)){
						z[i][j] = 4*p[i][j] - p[i][j-1] - p[i][j+1] - p[i-1][j];
					}else if (j==(n-1)){
						z[i][j] = 4*p[i][j] - p[i-1][j] - p[i+1][j] - p[i][j-1];
					}else{
						z[i][j] = 4*p[i][j] - p[i+1][j] - p[i-1][j] - p[i][j+1] - p[i][j-1];
					}
				}
			}
    		TAU_PROFILE_STOP(timer)			
			#pragma omp barrier

			#pragma omp single
			{
				vv=0.0;
			}
			#pragma omp barrier

    		TAU_PROFILE_TIMER(timer1, "dot_product1", "", TAU_USER);
    		TAU_PROFILE_START(timer1)
			#pragma omp for reduction(+: vv) collapse(2)  
			for(i = 0;i<n;i++){
				//#pragma omp parallel for  shared(i,n) num_threads(THREADS)
				for(j = 0;j<n;j++){
					vv += z[i][j]*p[i][j];
				}
			}
			TAU_PROFILE_STOP(timer1)			
			#pragma omp barrier
	    	#pragma omp single
	    	{
	    		alpha = del_old/vv;
	    	}
	    	#pragma omp barrier

    		TAU_PROFILE_TIMER(timer2, "vector_add1", "", TAU_USER);
    		TAU_PROFILE_START(timer2)
	    	#pragma omp for collapse(2)
			for(i = 0;i<n;i++){
				//#pragma omp parallel for schedule(static) shared(i,n) num_threads(THREADS)
				for(j = 0;j<n;j++){
					x[i][j] = x[i][j] + alpha*p[i][j];
				}
			}
			TAU_PROFILE_STOP(timer2)			
			#pragma omp barrier

    		TAU_PROFILE_TIMER(timer3, "vector_sub", "", TAU_USER);
    		TAU_PROFILE_START(timer3)
			#pragma omp for collapse(2) 
			for(i = 0;i<n;i++){
				//#pragma omp parallel for schedule(static) shared(i,n) num_threads(THREADS)
				for(j = 0;j<n;j++){
					r[i][j] = r[i][j] - alpha*z[i][j];
				}
			}
			TAU_PROFILE_STOP(timer3)			

			#pragma omp barrier
			
			#pragma omp single
			{
				del_new=0.0;
			}
			#pragma omp barrier
			
			TAU_PROFILE_TIMER(timer4, "dot_product2", "", TAU_USER);
    		TAU_PROFILE_START(timer4)
			#pragma omp for reduction(+: del_new) collapse(2)
			for(i = 0;i<n;i++){
				//#pragma omp parallel for  shared(i,n) num_threads(THREADS)
				for(j = 0;j<n;j++){
					del_new += r[i][j]*r[i][j];
				}
			}
			TAU_PROFILE_STOP(timer4)	

			#pragma omp barrier
			#pragma omp single
	    	{
	    		beta = del_new/del_old;
	    	}
	    	#pragma omp barrier

    		TAU_PROFILE_TIMER(timer5, "vector_add2", "", TAU_USER);
    		TAU_PROFILE_START(timer5)
	    	#pragma omp for collapse(2) 
			for(i = 0;i<n;i++){
				//#pragma omp parallel for schedule(static) shared(i,n) num_threads(THREADS)
				for(j = 0;j<n;j++){
					p[i][j] = r[i][j] + beta*p[i][j];
				}
			}
			TAU_PROFILE_STOP(timer5)			

			#pragma omp barrier

			#pragma omp single
	    	{
		    	del_old = del_new;
		    	float1 = del_new;
		    	float3 = sqrt(abs(float1));
				put = log10(float3/sqrt(float2) ) ;
		    	//cout << put << " ";
		    	iter++;
		    }
		    #pragma omp barrier

    	}
    	

    	
   }
    
    //cout << endl;
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>( t2 - t1 ).count();

    cout << "exec time:" << duration << endl;
    ofstream f;
    cout << "iterations: " << iter << endl;
    /*f.open("xtemp");
    rep(i,0,n-1){
    	rep(j,0,n-1){
    		f << std::fixed << std::setprecision(8) << x[i][j];
    		f << '\n';
    	}
    }
    f.close();*/
    free(b);
    free(r);
    free(z);
    free(x);
    free(p);
	

	

}