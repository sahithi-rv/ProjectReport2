#include<bits/stdc++.h>
#include <chrono>
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
#define output1(x) cout << #x << " " << x << endl;
#define output2(x,y) cout <<x<<" "<<y <<endl;
#define inp1(x) cin >> x;
#define inp2(x,y) cin >> x >> y;
#define MAX(a,b) a>b?a:b
#define MIN(a,b) a<b?a:b
#define N 1100


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

void matmul(float **b, float **out){
	rep(i,0,n-1){
		rep(j,0,n-1){
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
	rep(i,0,n-1){
		rep(j,0,n-1){
			out += l[i][j]*m[i][j];
		}
	}
	return out;
}

float dot_product1(float **l ){

	float out = 0.0;
	rep(i,0,n-1){
		rep(j,0,n-1){
			out += l[i][j]*l[i][j];
		}
	}
	return out;
}

float vector_add1(float **x, float alpha, float **p){
    rep(i,0,n-1){
		rep(j,0,n-1){
			x[i][j] = x[i][j] + alpha*p[i][j];
		}
	}
}

float vector_sub(float **r, float alpha, float **z){
	rep(i,0,n-1){
		rep(j,0,n-1){
			r[i][j] = r[i][j] - alpha*z[i][j];
		}
	}
}

float vector_add2(float **r, float beta, float **p){
   	rep(i,0,n-1){
		rep(j,0,n-1){
			p[i][j] = r[i][j] + beta*p[i][j];
		}
	}
}

int main(int argc, char *argv[]){

	//   	high_resolution_clock::time_point full1 = high_resolution_clock::now();
	if(argc<2){
		printf("Enter  arguments: input file for vector B,  Grid size\n");
		return 0;		
	}
	int iter = 0;
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
   			x[i][j] = 0.1*(rand()%10)+1;
   		}
   	}

   	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for(int i = 0;i<n;i++){
   		for(int j = 0;j<n;j++){
   			r[i][j] = b[i][j];
   			p[i][j] = b[i][j];
   			
        	del_new += b[i][j]*b[i][j];
        }
    }

    del_old = del_new;
	float float1 = del_new;
	float float2 = del_old;
	float float3 = sqrt(abs(float1));
	float put;
    //iterations
    float vv, alpha, beta;
    float leps = log10(EPS);
	put = log10(float3/sqrt(float2) ) ;
    while (put >= leps){
    	iter++;
    	matmul(p,z);

    	vv = dot_product2(z,p);
    	alpha = del_old/vv;

    	vector_add1(x,alpha, p);

    	vector_sub(r,alpha,z);


    	del_new = dot_product1(r);

    	beta = del_new/del_old;

    	vector_add2(r,beta,p);
 

    	del_old = del_new;
    	float1 = del_new;
    	float3 = sqrt(abs(float1));
		put = log10(float3/sqrt(float2) ) ;
    	//cout << put << " ";
   }
    
    //cout << endl;
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>( t2 - t1 ).count();
	 //  	high_resolution_clock::time_point full2 = high_resolution_clock::now();

	 //      auto duration2 = duration_cast<microseconds>( full2 - full1 ).count();

    cout << "exec time:" << duration  << endl;
    ofstream f;
    cout << "iterations: " << iter << endl;
    /*f.open("xtemp");
    rep(i,0,n-1){
    	rep(j,0,n-1){
    		f << std::fixed << std::setprecision(8) << x[i][j];
    		f << '\n';
    	}
    }
    f.close();
    free(b);
    free(r);
    free(z);
    free(x);
    free(p);*/

	

}