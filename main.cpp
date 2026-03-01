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
    double *a, *b, *x, *realx,*buf,*vecbuf,*resvec;
    // double el;
    double t1=0, t2=0;
    // feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);



    MPI_Init(&argc,&argv);

    int p;
    // MPI_Comm com;
    int main_kk = 0;
    int kk;
    // MPI_Comm_create()
    MPI_Comm_size(MPI_COMM_WORLD,&p);
    MPI_Comm_rank(MPI_COMM_WORLD,&kk);

    
    // MPI_Barrier(MPI_COMM_WORLD);

    filename=filename;

    if(!((argc == 5 || argc == 6) && sscanf(argv[1], "%d", &n)==1 && sscanf(argv[2], "%d", &m)==1  && sscanf(argv[3], "%d", &r)==1 && sscanf(argv[4], "%d", &s)==1)) 
    {
        if(kk == main_kk) cout<<"Usage : "<<argv[0]<<" <n> "<<" <m> " <<" <r> "<<" <s>"<<endl;
        MPI_Finalize();
        return 0;
    }
    // p = (n/m > p ? p:n/m);
   

    



    if(m<=0 || n<0 || r<0 || s<0 || p < 0)
    {
        if(kk == main_kk) printf("<n> or <m> or <r> or <s> or <p> <= 0, usage m,n,r,s,p > 0");
        MPI_Finalize();
        return 0;
    }

    if(argc == 6 && s!=0)
    {
        if(kk == main_kk) printf("Wrong usage! If s!=0 dont use initialization from file or s == 0 and file name not specified \n");
        // report(argv[0],task,r1,r2,t1,t2,s,n,m); 
        MPI_Finalize();
        return 0;
    }else if(argc == 5 && (s<1 || s>4))
    {
        if(kk == main_kk) printf("bad argument for s, s=1,2,3,4 \n");
        // report(argv[0],task,r1,r2,t1,t2,s,n,m); 
        MPI_Finalize();
        return 0;
    }

    //check
    if(argc == 6) 
    {
        filename = argv[5];
        // printf("n = %d m = %d r = %d s = %d p = %d k = %d file = %s\n",n,m,r,s,p,kk,filename);
    }
    // else 
    //     printf("n = %d m = %d r = %d s = %d p = %d k = %d\n",n,m,r,s,p,k);

    

   

    
    int loc_block_rows = get_block_rows(n,m,p,kk);
    loc_block_rows=loc_block_rows;
    int loc_rows = get_rows(n,m,p,kk);
    // printf("rows in %d proc %d\n",kk,loc_rows);
    
    // printf("a[%d,%d] = %10.3e\n",k,0,a[k*n+0]);
    a = new double[n*loc_rows]; //create matrix a
    b = new double[loc_rows];  // create vector b
    x = new double[loc_rows];  // create vector x
    realx = new double[loc_rows];  // create vector real x
    buf = new double[n*m]; //block row
    vecbuf = new double[m]; // block vector
    resvec = new double[m];

    // printf("b[0] = %10.3e proc %d\n",b[0],kk);
    
    for(int i_loc = 0; i_loc < loc_rows; i_loc++)
    {
        int i_glob = l2g(m,kk,p,i_loc);
        realx[i_loc] = (i_glob+1)%2;
        x[i_loc] = (i_glob+1)%2; // for template, later comment
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if(argc == 5) init_matrix(a,n,m,p,kk,&f,s);
    else if (argc == 6) 
    {
        if(read_matrix(a,n,m,p,kk,filename,buf,MPI_COMM_WORLD)!=0)
        {
            delete []a;
            delete []b;
            delete []x;
            delete []realx;
            delete []buf;
            delete []vecbuf;
            delete []resvec;
            r1 = -1; r2 = -1;
            if(kk == main_kk) report(argv[0],task,r1,r2,t1,t2,s,n,m,p);
            MPI_Finalize();
            return 0;
        }
    }

    if(kk == main_kk) printf("Matrix A:\n");
    print_matrix(a,n,m,p,kk,buf,r,MPI_COMM_WORLD);

    if(kk == main_kk) printf("Vector b:\n");
    MPI_Barrier(MPI_COMM_WORLD);
    init_vector_b(a,b,n,m,p,kk);
    print_vector(b,n,m,p,kk,vecbuf,r,MPI_COMM_WORLD);
    if(kk == main_kk) printf("\n\n");
    MPI_Barrier(MPI_COMM_WORLD);
    
    // print_array(b,n,1,1,r);

    // if(kk == main_kk) printf("Vector realx:\n");
    // print_vector(realx,n,m,p,kk,vecbuf,MPI_COMM_WORLD);

    // matrix_mult_vector(a,realx,resvec,n,m,kk,p,MPI_COMM_WORLD);

    // if(kk == main_kk) printf("Res mat_mult_vec:\n");
    // print_vector(resvec,n,m,p,kk,vecbuf,MPI_COMM_WORLD);
    
    int k = n/m;
    int l = n - k*m;
    
    double *block_mm = new double[m*m];
    double *block_ml = new double[m*l];
    double *block_ll = new double[l*l];
    double *tmpblock_mm = new double[m*m];
    double *tmpblock_ml = new double[m*l];
    double *tmpblock_ml1 = new double[m*l];
    double *tmpblock_ll = new double[l*l];
    double *invblock_mm = new double[m*m];
    double *invblock_ll = new double[l*l];
    double *diagblock_mm = new double[m*m];
    double *diaginvblock_mm = new double[m*m];
    double *vecb_m = new double[m];
    double *vecb_l = new double[l];
    double *tmpvecb_m = new double[m];
    double *tmpvecb_l = new double[l]; 
    /*static*/ int *colsw = new int[k];

    double elapsed = get_full_time();

    auto start_sol= std::chrono::high_resolution_clock::now();


    if(MPI_Solve(a,b,x,n,m,p,kk,buf,vecbuf,block_mm, block_ml, block_ll,invblock_mm, diaginvblock_mm, 
    invblock_ll,diagblock_mm,
    colsw,vecb_m,vecb_l,
    tmpblock_mm,tmpblock_ml,tmpblock_ml1,tmpblock_ll,tmpvecb_m,tmpvecb_l,MPI_COMM_WORLD) < 0)
    {
        delete []a;
        delete []b;
        delete []x;
        delete []realx;
        delete []buf;
        delete []vecbuf;
        delete []resvec;
        delete []block_mm ;
        delete []block_ml ;
        delete []block_ll ;
        delete []tmpblock_mm ;
        delete []tmpblock_ml ;
        delete []tmpblock_ml1 ;
        delete []tmpblock_ll ;
        delete []invblock_mm ;
        delete []invblock_ll ;
        delete []diagblock_mm ;
        delete []diaginvblock_mm ;
        delete []vecb_m ;
        delete []vecb_l ;
        delete []tmpvecb_m ;
        delete []tmpvecb_l ; 
        delete []colsw ;
        r1 = -1; r2 = -1;
        if(kk == main_kk) report(argv[0],task,r1,r2,t1,t2,s,n,m,p); // r1 = -1; r2 = -1;
        MPI_Finalize();
        return 0;
    }

    elapsed = get_full_time() - elapsed;
    printf("CPU Time proc %d = %.2lf\n",kk,elapsed);


    delete []block_mm ;
    delete []block_ml ;
    delete []block_ll ;
    delete []tmpblock_mm ;
    delete []tmpblock_ml ;
    delete []tmpblock_ml1 ;
    delete []tmpblock_ll ;
    delete []invblock_mm ;
    delete []invblock_ll ;
    delete []diagblock_mm ;
    delete []diaginvblock_mm ;
    delete []vecb_m ;
    delete []vecb_l ;
    delete []tmpvecb_m ;
    delete []tmpvecb_l ; 
    delete []colsw ;
    

    auto end_sol= std::chrono::high_resolution_clock::now();

    

    
    
    

    //print vector x
    
    
    

    auto start_res= std::chrono::high_resolution_clock::now();

    


    double *Ax = new double[loc_rows];//mat_x_vector(a,x,n);
    double *Ax_b = new double[loc_rows];//vectorsub( Ax , b, n);
    double *x_realx = new double[loc_rows];//vectorsub( x , realx, n);

    // for(int i = 0 ; i < n; i++)
    //     realx[i] = (i+1)%2;
    

    mpi_residuals(r1,r2,a,b,x,realx,n,m,p,kk,Ax,Ax_b,x_realx,MPI_COMM_WORLD);

    delete []Ax;
    delete []Ax_b;
    delete []x_realx;
    

    auto end_res= std::chrono::high_resolution_clock::now();

     

    t1 = chrono::duration<double>(end_sol - start_sol ).count();
    t2 = chrono::duration<double>(end_res - start_res).count();

    if(kk == main_kk) report(argv[0],task,r1,r2,t1,t2,s,n,m,p);

    



     

    delete []a;
    delete []b;
    delete []x;
    delete []realx;
    delete []buf;
    delete []vecbuf;
    delete []resvec;

    MPI_Finalize();


    return 0;
}
