#include <chrono>
#include <iostream>
#include "mpilib.h"




using namespace std;

int main(int argc, char *argv[])
{

    int n, m , r, s, task = 10;
    // int _c=0;
    char *filename = nullptr;
    double r1=0, r2=0;
    double *a, *b, *x, *realx;
    // double el;
    double t1=0, t2=0;
    // feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);



    MPI_Init(&argc,&argv);

    int p;
    // MPI_Comm com;
    int main_k = 0;
    int k;
    // MPI_Comm_create()
    MPI_Comm_size(MPI_COMM_WORLD,&p);
    MPI_Comm_rank(MPI_COMM_WORLD,&k);

    
    // MPI_Barrier(MPI_COMM_WORLD);

    filename=filename;

    if(!((argc == 5 || argc == 6) && sscanf(argv[1], "%d", &n)==1 && sscanf(argv[2], "%d", &m)==1  && sscanf(argv[3], "%d", &r)==1 && sscanf(argv[4], "%d", &s)==1)) 
    {
        cout<<"Usage : "<<argv[0]<<" <n> "<<" <m> " <<" <r> "<<" <s>"<<endl;
      
        return 0;
    }
    // p = (n/m > p ? p:n/m);
   

    



    if(m<=0 || n<0 || r<0 || s<0 || p < 0)
    {
        printf("<n> or <m> or <r> or <s> or <p> <= 0, usage m,n,r,s,p > 0");
        return 0;
    }

    if(argc == 6 && s!=0)
    {
        printf("Wrong usage! If s!=0 dont use initialization from file or s == 0 and file name not specified \n");
        // report(argv[0],task,r1,r2,t1,t2,s,n,m); 
        return 0;
    }else if(argc == 5 && (s<1 || s>4))
    {
        printf("bad argument for s, s=1,2,3,4 \n");
        // report(argv[0],task,r1,r2,t1,t2,s,n,m); 
        return 0;
    }

    //check
    // if(argc == 6) 
    // {
    //     filename = argv[5];
    //     printf("n = %d m = %d r = %d s = %d p = %d k = %d file = %s\n",n,m,r,s,p,k,filename);
    // }
    // else 
    //     printf("n = %d m = %d r = %d s = %d p = %d k = %d\n",n,m,r,s,p,k);

    

   

    a = new double[n*n]; //create matrix a
    b = new double[n];  // create vector b
    x = new double[n];  // create vector x
    realx = new double[n];  // create vector real x
    init_matrix(a,n,m,p,k,&f,s);
    // printf("a[%d,%d] = %10.3e\n",k,0,a[k*n+0]);
    
    
    

    // double elapsed = get_full_time();

    auto start_sol= std::chrono::high_resolution_clock::now();



    // printf("STATUS0 = %ld THREAD %d\n",status0,0);

     
    
    
        // printf("STATUS = %ld THREAD %d\n",status,thr);

       

        // elapsed = get_full_time() - elapsed;
        // printf("CPU Time thread %d = %.2lf\n",thr,elapsed);
    

    auto end_sol= std::chrono::high_resolution_clock::now();

    

    
    
    

    //print vector x
    
    
    

    auto start_res= std::chrono::high_resolution_clock::now();

    


    double *Ax = new double[n];//mat_x_vector(a,x,n);
    double *Ax_b = new double[n];//vectorsub( Ax , b, n);
    double *x_realx = new double[n];//vectorsub( x , realx, n);

    for(int i = 0 ; i < n; i++)
        realx[i] = (i+1)%2;
    

    // residuals(r1,r2,a,b,x,realx,n,Ax,Ax_b,x_realx);

    delete []Ax;
    delete []Ax_b;
    delete []x_realx;
    

    auto end_res= std::chrono::high_resolution_clock::now();

     

    t1 = chrono::duration<double>(end_sol - start_sol ).count();
    t2 = chrono::duration<double>(end_res - start_res).count();

    MPI_Finalize();



    report(argv[0],task,r1,r2,t1,t2,s,n,m,p); 

    delete []a;
    delete []b;
    delete []x;
    delete []realx;



    return 0;
}
