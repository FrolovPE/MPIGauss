#include <iostream>
#include <cmath>
#include "mpilib.h"
using namespace std;

void printlxn(const double *a, int size, int l, int n, int r)
{
    size=size;
    int l1 = min(l,r), n1 = min(n,r);


    for(int i = 0 ; i < l1 ; i++)
    {
        cout<<endl;
        for(int j =0  ;j < n1 ; j++)
        {
            printf("%10.3e ",a[i*n1+j]);
        }
    }
    cout<<endl;

}

int readarray(double *a, int n, const char* filename){
        int _c=0;
        double el;

        FILE *file = fopen(filename,"r");
        if(!file){
            
            printf("File %s doesnt exist or wrong file name!\n",filename);
            // report(argv[0],task,r1,r2,t1,t2,s,n,m); 
            // delete []a;
            // delete []b;
            // delete []x;
            // delete []realx;
            
            return -1;
        }

        while(fscanf(file,"%lf",&el)==1)
        {
            if(_c<n*n) 
            {
                a[_c] = el;
                _c++;
            }else
            {
                printf("Bad scan from file %s\n",filename);
                // report(argv[0],task,r1,r2,t1,t2,s,n,m); 
                // delete []a;
                // delete []b;
                // delete []x;
                // delete []realx;
                fclose(file);
                return -2;
            }

        
        }
        if(!feof(file) || _c!=n*n){
            printf("Bad file %s\n",filename);
            // report(argv[0],task,r1,r2,t1,t2,s,n,m); 
            fclose(file);
            // delete []a;
            // delete []b;
            // delete []x;
            // delete []realx;
            return -3;
        }
        fclose(file);
        return 0;
    }

int read_array(FILE *fp,double *buf,int &printed,int n,double totalsize)
{
    double x{};
    int counter{};
    if(!fp)
    {
        printf("cant open file\n");
        return -1;
    }

    while(counter < totalsize && fscanf(fp,"%lf",&x)==1)
    {
        buf[counter] = x;
        counter++;
        
    }
    printed += counter;
    if(printed == n*n && fscanf(fp,"%lf",&x)==1)
    {
        printf("bad size\n");
        return -1;
    }
    return 0;
}


double f (int s , int n , int i , int j)
{
    if(s == 1) 
        return n - max(i+1,j+1) +1;
    else if (s == 2)
        return max(i+1,j+1);
    else if (s == 3)
        return fabs(i-j);
    else if (s == 4)
        return static_cast<double>(1) / (static_cast<double>(i+1) + static_cast<double>(j+1) - static_cast<double>(1));
    else
        return -1;
}

void init(double *a,  double (*f)(int,int,int,int),  int n, int s )
{
    if(!a || !f)
    {
        printf("a or f = nullptr\n");
        return;
    }

    for(int i = 0; i < n ; i++)
    {
        for(int j=0 ; j < n; j++)
        {
            a[i*n+j] = f(s,n,i,j);
        }
    }
}

double vectornorm(double *a , int n)
{
    double res = 0;

    if(!a)
    {
        printf("nullptr in vector norm\n");
        return -1;
    }
    for(int i = 0; i<n; i++)
    {
        res += fabs(a[i]);
    }

    return res;
}

void mat_x_vector(double *res,double *a, double *b, int n)
{
    

    if(!a || !b)
    {
        return;
    }

    for(int i = 0; i< n; i++)
    {
        double s = 0;
        for(int j = 0; j < n; j++)
        {
            s += a[i*n +j] * b[j];
        }
        res[i] = s;
    }

  
}

void vectorsub(double *res,double *a,double *b, int n)
{
    // double *res = new double[n];

    if(!a || !b)
    {
        printf("nullptr in vector subtract\n");
        return;
    }

    for(int i = 0 ; i< n ; i++)
    {
        
        res[i] = a[i] - b[i];
    }

    // return res;
}



void residuals(double &r1,double &r2,double *a,double *b,double *x,double *realx,int n,double *Ax, double *Ax_b, double *x_realx)
{
    mat_x_vector(Ax,a,x,n);// double *Ax = mat_x_vector(a,x,n);
    vectorsub(Ax_b,Ax,b,n);// double *Ax_b = vectorsub( Ax , b, n);
    vectorsub(x_realx,x,realx,n);// double *x_realx = vectorsub( x , realx, n);

    // cout<<"\nVector Ax :"<<endl;
    // for(int i =0 ; i<n ; i++)
    // {
    //     printf("%10.3e ",Ax[i]);
    // }

    // cout<<"\nVector b :"<<endl;
    // for(int i =0 ; i<n ; i++)
    // {
    //     printf("%10.3e ",b[i]);
    // }

    // cout<<"\nVector Ax_b :"<<endl;
    // for(int i =0 ; i<n ; i++)
    // {
    //     printf("%10.3e ",Ax_b[i]);
    // }

    r1 = vectornorm(Ax_b ,  n) / vectornorm(b,n);
    r2 = vectornorm(x_realx ,  n) / vectornorm(realx,n);

    // delete []Ax;
    // delete []Ax_b;
    // delete []x_realx;
    
}

void mpi_residuals(double &r1,double &r2,double *a,double *b,double *x,double *realx,int n,int m,int p, int k, double *Ax, double *Ax_b, double *x_realx,MPI_Comm com)
{
    int lrows = get_rows(n,m,p,k);
    memset(Ax,0,lrows*sizeof(double));
    memset(Ax_b,0,lrows*sizeof(double));
    memset(x_realx,0,lrows*sizeof(double));
    /*dopisat' */matrix_mult_vector(a,x,Ax,n,m,k,p,com);// mat_x_vector(Ax,a,x,n);// double *Ax = mat_x_vector(a,x,n);
    // vectorsub(Ax_b,Ax,b,n);// double *Ax_b = vectorsub( Ax , b, n);
    // vectorsub(x_realx,x,realx,n);// double *x_realx = vectorsub( x , realx, n);
    double ln1{},ld1{},ln2{},ld2{},n1{},d1{},n2{},d2{};

    

    for(int i =0; i <lrows;i++)
    {
        ln1 += fabs(Ax[i]-b[i]);
        ld1+=fabs(b[i]);

        ln2+=fabs(x[i]-realx[i]);
        ld2+=(realx[i]);
    }

    MPI_Allreduce(&ln1,&n1,1,MPI_DOUBLE,MPI_SUM,com);
    MPI_Allreduce(&ln2,&n2,1,MPI_DOUBLE,MPI_SUM,com);
    MPI_Allreduce(&ld1,&d1,1,MPI_DOUBLE,MPI_SUM,com);
    MPI_Allreduce(&ld2,&d1,1,MPI_DOUBLE,MPI_SUM,com);

    r1 = (d1 > 0) ? n1/d1:0;
    r2 = (d2>0) ? n2/d2:0;

    // cout<<"\nVector Ax :"<<endl;
    // for(int i =0 ; i<n ; i++)
    // {
    //     printf("%10.3e ",Ax[i]);
    // }

    // cout<<"\nVector b :"<<endl;
    // for(int i =0 ; i<n ; i++)
    // {
    //     printf("%10.3e ",b[i]);
    // }

    // cout<<"\nVector Ax_b :"<<endl;
    // for(int i =0 ; i<n ; i++)
    // {
    //     printf("%10.3e ",Ax_b[i]);
    // }



    // delete []Ax;
    // delete []Ax_b;
    // delete []x_realx;
    
}






void report(char *title, int task, double r1, double r2 ,double t1,  double t2 ,int s, int n , int m, int p )
{
    printf ("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2lf T2 = %.2lf S = %d N = %d M = %d P = %d\n",
title, task, r1, r2, t1, t2, s, n, m,p);
}


void matmult(double *res,double *a, double *b, int n, int m,int l) //a^{n*m} * b^{m*l} = res^{n*l}
{
    int i = 0,j=0;
    
   
    if(!a || !b)
    {
        return;
    }

    for(i = 0 ; i < n ; i++)
    {
        for(j = 0; j< l; j++)
        {
            double s = 0;
            for(int k = 0 ; k < m ; k++)
            {
                s += a[i*m+k]*b[k*l+j];
            }
            res[i*l+j] = s;
        }
    }

}

 void get_block(double *a, double *b, int n, int m, int i, int j)
{
    int i1=0, j1=0, k, l, r, h;
    if(m == 0) {
        printf("m == 0 in i = %d , j = %d\n",i,j);
        return;
    }
    k = n/m; l = n - m*k ;
    if(i < k) r = m;
    else r = l;

    if(j < k) h = m;
    else h = l;


    double *bl = a + i*n*m + j*m; //start of block double *bl = a + i*n*m + j*m;

    for (i1 =0 ; i1 < r ; i1++)
    {
        for (j1=0 ; j1 < h ; j1++)
        {
            b[i1*h + j1] = bl[i1*n + j1];
            

        }
    }
    

}



void set_block(double *a, double *b, int n, int m, int i, int j)
{
    int i1=0, j1=0, k, l, r, h;
    k = n/m; l = n - m*k ;
    if(i < k) r = m;
    else r = l;

    if(j < k) h = m;
    else h = l;


    double *bl = a + i*n*m + j*m; //start of block

    for (i1 =0 ; i1 < r ; i1++)
    {
        for (j1=0 ; j1 < h ; j1++)
        {
            bl[i1*n + j1] = b[i1*h + j1];
            

        }
    }
    

}

double normofmatrix(double *a , int size)
{
    if (!a)
    { 
        printf("nullptr in norm of matrix\n");
        return -1;
    }

    double mm = -1;

    for(int j = 0; j< size ; j++)
    {
        double s = 0;
        for(int i = 0 ;i <size; i++)
        {
            s+=fabs(a[i*size+j]);
        }
        if(mm < s) mm = s;
    }

    return mm;

}

void mpi_norm_matrix(double *a,double &N,int n,int m, int p, int k,MPI_Comm com)
{
    double norm = -1;
    
    int rows = get_rows(n,m,p,k);
    int i_loc,j_loc;

    for(i_loc = 0; i_loc < rows; i_loc++)
    {double sum = 0;
        for(j_loc = 0; j_loc < n; j_loc++)
        {
            sum += fabs(a[i_loc*n + j_loc]);
        }
        if(norm < sum) 
        {
            norm = sum;
            // printf("proc %d norm = %10.3e\n",k,norm);
        }
    }

    MPI_Allreduce(&norm,&N,1,MPI_DOUBLE,MPI_MAX,com);
}

double* inverse(double *result,const double* A, int size,double eps)
{
    if(!A)
    {
        printf("nullptr in inverse\n");
        return nullptr;
    }

    // for(int i = 0; i<size*size; i++)
    // {
    //     if(A[i])
    //     {
    //         printf("cant inverse non square matrix\n");
    //         return nullptr;
    //     }
    // }

    double* a = new double[size*size];
    memcpy(a, A, size*size * sizeof(double));


    double* E = new double[size*size];
    int* colsw = new int[size];
    
    // init E
    for(int i = 0; i<size; i++)
    {
        for(int j = 0; j< size; j++)
        {
            if (i!=j) E[i*size+j] = 0;
            else E[i*size +j] = 1;
        }
    }

    //init colsw
    for(int i = 0; i <size ;i++) colsw[i] = i;


    for(int i = 0 ; i< size ; i++)
    {
        double maxEl = fabs(a[i*size+i]);
        int pivot = i;
        for(int j = i+1; j<size;j++)
        {
            if(fabs(a[i*size +j]) > maxEl)
            {
                maxEl = fabs(a[i*size +j]);
                pivot = j;
            }
        }

        
        //swap columns
        if(pivot!=i)
        {
            for(int k = 0; k<size ; k++)
            {
                swap(a[k*size+i],a[k*size+pivot]);
                swap(E[k*size+i],E[k*size+pivot]);
            }
            swap(colsw[i],colsw[pivot]);
        }

        if(fabs(a[i*size+i]) < eps)
        {
            // printf("matrix has no inverse\n");
            delete []E;
            delete []colsw;
            delete []a;
            return nullptr;
        }

        //devide row i 
        // cout<<"A matr before devideing"<<endl;
        // printlxn(a,size,size,size,size);

        double mainEl = a[i*size+i];
        for(int k = 0 ; k<size; k++)
        {
            a[i*size+k] /= mainEl;
            E[i*size+k] /= mainEl;
        }
        // cout<<"A matr after devideing"<<endl;
        // printlxn(a,size,size,size,size);


        for(int k = 0; k< size;k++)
        {
            if(k!=i)
            {
                double factor = a[k*size+i];
                for(int j = 0; j <size;j++)
                {
                    a[k*size+j] -= factor * a[i*size+j];
                    E[k*size+j] -= factor * E[i*size+j];
                }
            }
        }

    }

// cout<<"Last a presentation before swap"<<endl;
// printlxn(a,size,size,size,size);


    // swap back columns
    

    for (int i = 0; i < size*size; ++i) result[i] = 0.0;

    
    
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            result[ colsw[i]*size + colsw[j] ] = E[i*size + j];
// cout<<"Last a presentation after swap"<<endl;

// printlxn(a,size,size,size,size);

    delete []a;
    delete [] colsw;
    delete []E;
    return result;
    
    
}

void blocksize(int i, int j,int n,int m, int &r, int &h)
{
    int k,l;
    k = n/m; l = n - m*k ;
    if(i < k) r = m;
        else r = l;

    if(j < k) h = m;
    else h = l;
}


void swap_block_columns(double *a, int n,int m, int i, int j)
{
    for(int p =0 ; p < m ; p++)
                {
                    for(int c = 0; c<n ; c++)
                    {
                        //должны свапнуть столбцы блоков 
                        swap(a[c*n+i*m+p],a[c*n+j*m+p]);
                        
                    }
                }
                
                // printf("swapped blocks %d and %d\n",i,j);
}

void swap_block_vec(double *a, int n,int m, int i, int j)
{
    n=n;
    for(int p =0 ; p < m ; p++)
                {
                    
                        //должны свапнуть столбцы блоков 
                        swap(a[i*m+p],a[j*m+p]);
                        
                    
                }
                
               
}

void get_vec_block(double *b,double *block,int n, int m, int i)
{
    int i1=0, k, l, r;
    k = n/m; l = n - m*k ;
    if(i < k) r = m;
    else r = l;

    for (i1 =0 ; i1 < r ; i1++)
    {
            block[i1] = b[i*m+i1];
    }
    

}


void mat_mult_sub(double* Result, double* Block_A, double* Block_B, const int row_A, const int col_B, const int col_row) {
    int row_l = row_A % 3;
    int col_l = col_B % 3;
    int row_k = (row_A - row_l) / 3;
    int col_k = (col_B - col_l) / 3;
    double res_00 = 0, res_01 = 0, res_02 = 0;
    double res_10 = 0, res_11 = 0, res_12 = 0;
    double res_20 = 0, res_21 = 0, res_22 = 0;

    for (int b_i = 0; b_i < row_k; b_i++) {
        for (int b_j = 0; b_j < col_k; b_j++) {
            res_00 = 0, res_01 = 0, res_02 = 0;
            res_10 = 0, res_11 = 0, res_12 = 0;
            res_20 = 0, res_21 = 0, res_22 = 0;
            for (int s = 0; s < col_row; s++) {
                res_00 += Block_A[b_i * 3 * col_row + s] * Block_B[s * col_B + b_j * 3];
                res_01 += Block_A[b_i * 3 * col_row + s] * Block_B[s * col_B + b_j * 3 + 1];
                res_02 += Block_A[b_i * 3 * col_row + s] * Block_B[s * col_B + b_j * 3 + 2];

                res_10 += Block_A[(b_i * 3 + 1) * col_row + s] * Block_B[s * col_B + b_j * 3];
                res_11 += Block_A[(b_i * 3 + 1) * col_row + s] * Block_B[s * col_B + b_j * 3 + 1];
                res_12 += Block_A[(b_i * 3 + 1) * col_row + s] * Block_B[s * col_B + b_j * 3 + 2]; 

                res_20 += Block_A[(b_i * 3 + 2) * col_row + s] * Block_B[s * col_B + b_j * 3];
                res_21 += Block_A[(b_i * 3 + 2) * col_row + s] * Block_B[s * col_B + b_j * 3 + 1];
                res_22 += Block_A[(b_i * 3 + 2) * col_row + s] * Block_B[s * col_B + b_j * 3 + 2];
            }

            Result[b_i * 3 * col_B + b_j * 3] -= res_00; Result[b_i * 3 * col_B + b_j * 3 + 1] -= res_01; Result[b_i * 3 * col_B + b_j * 3 + 2] -= res_02;
            Result[(b_i * 3 + 1) * col_B + b_j * 3] -= res_10; Result[(b_i * 3 + 1) * col_B + b_j * 3 + 1] -= res_11; Result[(b_i * 3 + 1) * col_B + b_j * 3 + 2] -= res_12;
            Result[(b_i * 3 + 2) * col_B + b_j * 3] -= res_20; Result[(b_i * 3 + 2) * col_B + b_j * 3 + 1] -= res_21; Result[(b_i * 3 + 2) * col_B + b_j * 3 + 2] -= res_22;
        }

        if (col_l != 0) {
            res_00 = 0, res_01 = 0, res_10 = 0;
            res_11 = 0, res_20 = 0, res_21 = 0;

            for (int s = 0; s < col_row; s++) {
                if(col_l > 1) {
                    res_01 += Block_A[b_i * 3 * col_row + s] * Block_B[s * col_B + col_k * 3 + 1];
                    res_11 += Block_A[(b_i * 3 + 1) * col_row + s] * Block_B[s * col_B + col_k * 3 + 1];
                    res_21 += Block_A[(b_i * 3 + 2) * col_row + s] * Block_B[s * col_B + col_k * 3 + 1];
                }

                res_00 += Block_A[b_i * 3 * col_row + s] * Block_B[s * col_B + col_k * 3];
                res_10 += Block_A[(b_i * 3 + 1) * col_row + s] * Block_B[s * col_B + col_k * 3];
                res_20 += Block_A[(b_i * 3 + 2) * col_row + s] * Block_B[s * col_B + col_k * 3];
            }

            Result[b_i * 3 * col_B + col_k * 3] -= res_00; 
            Result[(b_i * 3 + 1) * col_B + col_k * 3] -= res_10;
            Result[(b_i * 3 + 2) * col_B + col_k * 3] -= res_20;

            if(col_l > 1) {
                Result[b_i * 3 * col_B + col_k * 3 + 1] -= res_01; 
                Result[(b_i * 3 + 1) * col_B + col_k * 3 + 1] -= res_11;
                Result[(b_i * 3 + 2) * col_B + col_k * 3 + 1] -= res_21;
            }
        }
            
    }

    if(row_l != 0) {
        for (int b_j = 0; b_j < col_k; b_j++) {
            res_00 = 0, res_01 = 0, res_02 = 0;
            res_10 = 0, res_11 = 0, res_12 = 0;
            for (int s = 0; s < col_row; s++) {
                res_00 += Block_A[(row_k * 3) * col_row + s] * Block_B[s * col_B + b_j * 3];
                res_01 += Block_A[(row_k * 3) * col_row + s] * Block_B[s * col_B + b_j * 3 + 1];
                res_02 += Block_A[(row_k * 3) * col_row + s] * Block_B[s * col_B + b_j * 3 + 2]; 

                if(row_l > 1) {
                    res_10 += Block_A[(row_k * 3 + 1) * col_row + s] * Block_B[s * col_B + b_j * 3];
                    res_11 += Block_A[(row_k * 3 + 1) * col_row + s] * Block_B[s * col_B + b_j * 3 + 1];
                    res_12 += Block_A[(row_k * 3 + 1) * col_row + s] * Block_B[s * col_B + b_j * 3 + 2]; 
                }
            }

            Result[row_k * 3 * col_B + b_j * 3] -= res_00; Result[row_k * 3 * col_B + b_j * 3 + 1] -= res_01; Result[row_k * 3 * col_B + b_j * 3 + 2] -= res_02;

            if (row_l > 1) {
                Result[(row_k * 3 + 1) * col_B + b_j * 3] -= res_10; Result[(row_k * 3 + 1) * col_B + b_j * 3 + 1] -= res_11; Result[(row_k * 3 + 1) * col_B + b_j * 3 + 2] -= res_12;
            }
        }

        if(col_l != 0) {
            res_00 = 0, res_01 = 0;
            res_10 = 0, res_11 = 0;
            for (int s = 0; s < col_row; s++) {
                res_00 += Block_A[row_k * 3 * col_row + s] * Block_B[s * col_B + col_k * 3];
                if (col_l > 1) {
                    res_01 += Block_A[row_k * 3 * col_row + s] * Block_B[s * col_B + col_k * 3 + 1];
                }
                if(row_l > 1) {
                    res_10 += Block_A[(row_k * 3  + 1) * col_row + s] * Block_B[s * col_B + col_k * 3];
                }
                if (col_l > 1 && row_l > 1) {
                    res_11 += Block_A[(row_k * 3  + 1) * col_row + s] * Block_B[s * col_B + col_k * 3 + 1];
                }
            }

            Result[row_k * 3 * col_B + col_k * 3] -= res_00;

            if (col_l > 1) {
                Result[row_k * 3 * col_B + col_k * 3 + 1] -= res_01;
            }
            
            if (row_l > 1) {
                Result[(row_k * 3 + 1) * col_B + col_k * 3] -= res_10;
            }

            if (row_l > 1 && col_l > 1) {
                Result[(row_k * 3 + 1) * col_B + col_k * 3 + 1] -= res_11;
            }
        }     
    }
}

void set_vec_block(double *b,double *block,int n, int m, int i)
{
    int i1=0, k, l, r;
    k = n/m; l = n - m*k ;
    if(i < k) r = m;
    else r = l;

    


    //start of block 

    for (i1 =0 ; i1 < r ; i1++)
    {
            b[i*m+i1]=block[i1] ;
    }
    

}

void get_block_ml(double *a, double *b, int n, int m,int l, int i)
{
    
    int i1=0, j1=0;
   

    double *bl = a + i*n*m + (n-l); //start of block double *bl = a + i*n*m + j*m;

    for (i1 =0 ; i1 < m ; i1++)
    {
        for (j1=0 ; j1 < l ; j1++)
        {
            b[i1*l + j1] = bl[i1*n + j1];
            

        }
    }
    

    
}

void set_block_ml(double *a, double *b, int n, int m,int l, int i)
{
    
    int i1=0, j1=0;
   

    double *bl = a + i*n*m + (n-l); //start of block double *bl = a + i*n*m + j*m;

    for (i1 =0 ; i1 < m ; i1++)
    {
        for (j1=0 ; j1 < l ; j1++)
        {
             bl[i1*n + j1]=b[i1*l + j1] ;
            

        }
    }
    

    
}

void get_block_lm(double *a, double *b, int n, int m,int l, int j)
{
int i1=0, j1=0;
   

    double *bl = a + (n-l)*n + j*m; //start of block double *bl = a + i*n*m + j*m;

    for (i1 =0 ; i1 < l ; i1++)
    {
        for (j1=0 ; j1 < m ; j1++)
        {
            b[i1*m + j1] = bl[i1*n + j1];
            

        }
    }
    
}

void set_block_lm(double *a, double *b, int n, int m,int l, int j)
{
int i1=0, j1=0;
   

    double *bl = a + (n-l)*n + j*m; //start of block double *bl = a + i*n*m + j*m;

    for (i1 =0 ; i1 < l ; i1++)
    {
        for (j1=0 ; j1 < m ; j1++)
        {
            bl[i1*n + j1]= b[i1*m + j1];
            

        }
    }
    
}

void matsub(double *res,double *a, double *b, int m,int l)
{
    for(int i = 0 ;i<m;i++)
    {
        for(int j = 0; j<l; j++)
        {
            res[i*l+j] = a[i*l+j] - b[i*l+j];
        }
    }
}

void vec_mult_sub(double* Result, double* A, double* vec, int m) {
    double* temp = new double[m];

    for (int i = 0; i < m; i++)
    {
        
        temp[i] = 0;
        double sum = 0;
        for (int j = 0; j < m; j++) {
            
            if(fabs( A[i*m + j] ) < EPS64 ) 
            {
                A[i*m + j] = 0;
            }

                sum += A[i*m + j]*vec[j];

        }
        temp[i] += sum;
            // printf("temp[%d] = %lf\n",i,temp[i]);
        
    }

    for (int i = 0; i < m; i++)
        Result[i] -= temp[i];

    // cout<<"vec_mult:\n";
    // for(int i = 0 ; i < m ; i++)
    // {
    //     cout<<Result[i]<<" ";
    // }
    delete[] temp;
}

void vec_mult_sub_lm(double* Result, double* A, double* vec, int l,int m) {
    double* temp = new double[l];
    for (int i = 0; i < l; i++)
    {
        temp[i] = 0;
        double sum = 0;
        for (int j = 0; j < m; j++) {
            
            if(fabs(A[i*m + j]) < EPS64 ) 
            {
                A[i*m + j] = 0;
            }
                sum += A[i*m + j]*vec[j];

        }
        temp[i] += sum;
            // printf("temp[%d] = %lf\n",i,temp[i]);
        
    }

    for (int i = 0; i < l; i++)
        Result[i] -= temp[i];

    // cout<<"vec_mult:\n";
    // for(int i = 0 ; i < m ; i++)
    // {
    //     cout<<Result[i]<<" ";
    // }
    delete[] temp;
}







void multiplication(double* Result, double* Block_A, double* Block_B, const int row_A, const int col_row,const int col_B)
{
    int row_l = row_A % 3;
    int col_l = col_B % 3;
    int row_k = (row_A - row_l) / 3;
    int col_k = (col_B - col_l) / 3;
    double res_00 = 0, res_01 = 0, res_02 = 0;
    double res_10 = 0, res_11 = 0, res_12 = 0;
    double res_20 = 0, res_21 = 0, res_22 = 0;

    for (int b_i = 0; b_i < row_k; b_i++) {
        for (int b_j = 0; b_j < col_k; b_j++) {
            res_00 = 0, res_01 = 0, res_02 = 0;
            res_10 = 0, res_11 = 0, res_12 = 0;
            res_20 = 0, res_21 = 0, res_22 = 0;
            for (int s = 0; s < col_row; s++) {
                res_00 += Block_A[b_i * 3 * col_row + s] * Block_B[s * col_B + b_j * 3];
                res_01 += Block_A[b_i * 3 * col_row + s] * Block_B[s * col_B + b_j * 3 + 1];
                res_02 += Block_A[b_i * 3 * col_row + s] * Block_B[s * col_B + b_j * 3 + 2];

                res_10 += Block_A[(b_i * 3 + 1) * col_row + s] * Block_B[s * col_B + b_j * 3];
                res_11 += Block_A[(b_i * 3 + 1) * col_row + s] * Block_B[s * col_B + b_j * 3 + 1];
                res_12 += Block_A[(b_i * 3 + 1) * col_row + s] * Block_B[s * col_B + b_j * 3 + 2]; 

                res_20 += Block_A[(b_i * 3 + 2) * col_row + s] * Block_B[s * col_B + b_j * 3];
                res_21 += Block_A[(b_i * 3 + 2) * col_row + s] * Block_B[s * col_B + b_j * 3 + 1];
                res_22 += Block_A[(b_i * 3 + 2) * col_row + s] * Block_B[s * col_B + b_j * 3 + 2];
            }

            Result[b_i * 3 * col_B + b_j * 3] = res_00; Result[b_i * 3 * col_B + b_j * 3 + 1] = res_01; Result[b_i * 3 * col_B + b_j * 3 + 2] = res_02;
            Result[(b_i * 3 + 1) * col_B + b_j * 3] = res_10; Result[(b_i * 3 + 1) * col_B + b_j * 3 + 1] = res_11; Result[(b_i * 3 + 1) * col_B + b_j * 3 + 2] = res_12;
            Result[(b_i * 3 + 2) * col_B + b_j * 3] = res_20; Result[(b_i * 3 + 2) * col_B + b_j * 3 + 1] = res_21; Result[(b_i * 3 + 2) * col_B + b_j * 3 + 2] = res_22;
        }

        if (col_l != 0) {
            res_00 = 0, res_01 = 0, res_10 = 0;
            res_11 = 0, res_20 = 0, res_21 = 0;

            for (int s = 0; s < col_row; s++) {
                if(col_l > 1) {
                    res_01 += Block_A[b_i * 3 * col_row + s] * Block_B[s * col_B + col_k * 3 + 1];
                    res_11 += Block_A[(b_i * 3 + 1) * col_row + s] * Block_B[s * col_B + col_k * 3 + 1];
                    res_21 += Block_A[(b_i * 3 + 2) * col_row + s] * Block_B[s * col_B + col_k * 3 + 1];
                }

                res_00 += Block_A[b_i * 3 * col_row + s] * Block_B[s * col_B + col_k * 3];
                res_10 += Block_A[(b_i * 3 + 1) * col_row + s] * Block_B[s * col_B + col_k * 3];
                res_20 += Block_A[(b_i * 3 + 2) * col_row + s] * Block_B[s * col_B + col_k * 3];
            }

            Result[b_i * 3 * col_B + col_k * 3] = res_00; 
            Result[(b_i * 3 + 1) * col_B + col_k * 3] = res_10;
            Result[(b_i * 3 + 2) * col_B + col_k * 3] = res_20;

            if(col_l > 1) {
                Result[b_i * 3 * col_B + col_k * 3 + 1] = res_01; 
                Result[(b_i * 3 + 1) * col_B + col_k * 3 + 1] = res_11;
                Result[(b_i * 3 + 2) * col_B + col_k * 3 + 1] = res_21;
            }
        }
            
    }

    if(row_l != 0) {
        for (int b_j = 0; b_j < col_k; b_j++) {
            res_00 = 0, res_01 = 0, res_02 = 0;
            res_10 = 0, res_11 = 0, res_12 = 0;
            for (int s = 0; s < col_row; s++) {
                res_00 += Block_A[(row_k * 3) * col_row + s] * Block_B[s * col_B + b_j * 3];
                res_01 += Block_A[(row_k * 3) * col_row + s] * Block_B[s * col_B + b_j * 3 + 1];
                res_02 += Block_A[(row_k * 3) * col_row + s] * Block_B[s * col_B + b_j * 3 + 2]; 

                if(row_l > 1) {
                    res_10 += Block_A[(row_k * 3 + 1) * col_row + s] * Block_B[s * col_B + b_j * 3];
                    res_11 += Block_A[(row_k * 3 + 1) * col_row + s] * Block_B[s * col_B + b_j * 3 + 1];
                    res_12 += Block_A[(row_k * 3 + 1) * col_row + s] * Block_B[s * col_B + b_j * 3 + 2]; 
                }
            }

            Result[row_k * 3 * col_B + b_j * 3] = res_00; Result[row_k * 3 * col_B + b_j * 3 + 1] = res_01; Result[row_k * 3 * col_B + b_j * 3 + 2] = res_02;

            if (row_l > 1) {
                Result[(row_k * 3 + 1) * col_B + b_j * 3] = res_10; Result[(row_k * 3 + 1) * col_B + b_j * 3 + 1] = res_11; Result[(row_k * 3 + 1) * col_B + b_j * 3 + 2] = res_12;
            }
        }

        if(col_l != 0) {
            res_00 = 0, res_01 = 0;
            res_10 = 0, res_11 = 0;
            for (int s = 0; s < col_row; s++) {
                res_00 += Block_A[row_k * 3 * col_row + s] * Block_B[s * col_B + col_k * 3];
                if (col_l > 1) {
                    res_01 += Block_A[row_k * 3 * col_row + s] * Block_B[s * col_B + col_k * 3 + 1];
                }
                if(row_l > 1) {
                    res_10 += Block_A[(row_k * 3  + 1) * col_row + s] * Block_B[s * col_B + col_k * 3];
                }
                if (col_l > 1 && row_l > 1) {
                    res_11 += Block_A[(row_k * 3  + 1) * col_row + s] * Block_B[s * col_B + col_k * 3 + 1];
                }
            }

            Result[row_k * 3 * col_B + col_k * 3] = res_00;

            if (col_l > 1) {
                Result[row_k * 3 * col_B + col_k * 3 + 1] = res_01;
            }
            
            if (row_l > 1) {
                Result[(row_k * 3 + 1) * col_B + col_k * 3] = res_10;
            }

            if (row_l > 1 && col_l > 1) {
                Result[(row_k * 3 + 1) * col_B + col_k * 3 + 1] = res_11;
            }
        }     
    }      
}



int solution(int n, int m, double *a, double *b, double *x,
    double *block_mm, double *block_ml, double *block_ll,double *invblock_mm, double *diaginvblock_mm, 
    double *invblock_ll,double *diagblock_mm,
    int *colsw,double *vecb_m,double *vecb_l,
    double *tmpblock_mm,double *tmpblock_ml,double *tmpblock_ml1,double *tmpblock_ll,double *tmpvecb_m,double *tmpvecb_l,const double eps)
{
    
    // m=m;
    // a=a;
    // b=b;
    // block_ml=block_ml;
    // diaginvblock_mm=diaginvblock_mm;
    // diagblock_mm=diagblock_mm;
    // vecb_l=vecb_l;
    // tmpblock_ml1 = tmpblock_ml1;

    

    

    int  k, l/*, r, h*/;
    // int i = 1,
    // j=0;

    k = n/m; l = n - m*k ;
    
    int is_l = (l == 0) ? 0:1; 
    
    
    for(int c =0; c < k  ;c++) colsw[c]=c;

    for(int i = 0 ; i < k + is_l; i++)
    {   
        double minNorm = 1e64;

        int mainBlock = i;

        if(i != k)
        {
            for(int j = i ; j< k ; j++)
            {

            get_block(a,block_mm,n,m,i,j);

    //             printf("Block[%d,%d]\n",i,j);
    //             printlxn(block_mm,m,m,m,m);

            
            // printlxn(invblock_mm,m,m,m,m);

            if(inverse(invblock_mm,block_mm,m,eps))
                {

                        
//                     cout<<"inverse "<<i<<" "<<j<<" with norm = "<<normofmatrix(invblock_mm,m)<<endl;
// 
//                 printlxn(invblock_mm,m,m,m,m);

                if(normofmatrix(invblock_mm,m) < minNorm) 
                    {
                        
                        minNorm = normofmatrix(invblock_mm,m);
                        mainBlock = j;
                       
                    }
                }

            }

        }else{
            get_block(a,block_ll,n,m,k,k);

            if(!inverse(invblock_ll,block_ll,l,eps))
            {
                printf("Block [%d,%d] (block[l,l] in our matrix)  has no inverse after the transformations\n",k,k);
                return -1;
            }
            
            
//             printlxn(invblock_ll,l,l,l,l);
            minNorm = normofmatrix(invblock_ll,l);
            
        }

        if((fabs(minNorm - 1e64) < eps))
        {   
            if(i!=0)
                printf("No inverse matrix in row %d after the transformations\n",i);
            else
                printf("No inverse matrix in row %d\n",i);
            return -1;
        }

        if(mainBlock != i)
            {
                swap_block_columns(a,n,m,i,mainBlock);
                // printlxn(a,n,n,n,n);
                swap(colsw[i],colsw[mainBlock]);
                // cout<<"swapped "<< i<<" "<<mainBlock<<" in row "<<i<<endl;
            }
        
        
        // printlxn(a,n,n,n,n);
        // cout<<"TEST1"<<endl;
        if(i<k)
        {
            get_block(a,diagblock_mm,n,m,i,i);
            
            if(!(inverse(diaginvblock_mm,diagblock_mm,m,eps)))
            {
                        printf("no blocks in row has inverse block\n");
                        
                        return -1;
                    }

            get_vec_block(b,vecb_m,n,m,i);
            mat_x_vector(tmpvecb_m,diaginvblock_mm,vecb_m,m);// double *resvec = mat_x_vector(diaginvblock_mm,vecb_m,m);
            // cout<<"tmpvecb_m : "<<endl;
            // printlxn(tmpvecb_m,m,1,m,m);    
            set_vec_block(b,tmpvecb_m,n,m,i);

            for(int j = i ; j < k ; j++) //mb try j = i
            {
                get_block(a,block_mm,n,m,i,j);
                
               multiplication(tmpblock_mm,diaginvblock_mm,block_mm,m,m,m);// matmult(tmpblock_mm,diaginvblock_mm,block_mm,m,m,m);// double *resmult = matmult(diaginvblock_mm,block_mm,m,m,m)

                set_block(a,tmpblock_mm,n,m,i,j);
                
                            if (!block_mm || !vecb_m || !invblock_mm) {
                fprintf(stderr, "Error: temporary buffers not initialized!\n");
                return -1;
            }
            }
            // printlxn(a,n,n,n,n);
            if(is_l != 0)
            {
                get_block_ml(a,block_ml,n,m,l,i);
                multiplication(tmpblock_ml,diaginvblock_mm,block_ml,m,m,l);// matmult(tmpblock_ml,diaginvblock_mm,block_ml,m,m,l);
                set_block_ml(a,tmpblock_ml,n,m,l,i);
            }
            
        }else
            {   
                // printlxn(a,n,n,n,n);
                // printlxn(b,n,1,n,n);
                get_block(a,block_ll,n,m,i,i);
                get_vec_block(b,vecb_l,n,m,i);

                // printlxn(block_ll,l,l,l,n);

                // cout<<"vecb_l:"<<endl;
                // printlxn(vecb_l,l,1,l,n);

                if(!(inverse(invblock_ll,block_ll,l,eps)))
                    {
                        printf("ll block has no inverse\n");
                        return -1;
                    }

                
                // printlxn(invblock_ll,l,l,l,n);

                multiplication(tmpblock_ll,invblock_ll,block_ll,l,l,l);// matmult(tmpblock_ll,invblock_ll,block_ll,l,l,l);

                mat_x_vector(tmpvecb_l,invblock_ll,vecb_l,l);
                
                // printlxn(tmpvecb_l,l,1,l,n);

                set_block(a,tmpblock_ll,n,m,i,i);
                set_vec_block(b,tmpvecb_l,n,m,i);
                // cout<<"WE ARE IN i = k"<<endl;
            }

            // printlxn(a,n,n,n,n);
            // printlxn(b,n,1,n,n);
            //начинаем обнулять столбцы
            
        for(int r = i+1 ; r < k + is_l ; r++)
        {
            // cout<<"TEST "<<r<<endl;
            if(r < k)
            {
                get_block(a,block_mm,n,m,r,i);
                get_block(a,tmpblock_mm,n,m,r,i);
                memset(tmpblock_mm,0, m*m*sizeof(double));
                set_block(a,tmpblock_mm,n,m,r,i);

                // not in i for
                get_vec_block(b,vecb_m,n,m,i);//вычитание из вектора b block_mm*b
                get_vec_block(b,tmpvecb_m,n,m,r);
                vec_mult_sub(tmpvecb_m,block_mm,vecb_m,m);
                set_vec_block(b,tmpvecb_m,n,m,r);

                // cout<<"tmpvecb_m in subtract i= "<<i<<" r="<<r<<endl;
                // printlxn(tmpvecb_m,m,1,m,m);
                // printlxn(b,n,1,n,n);

                for (int j = i + 1; j < k; j++) {
                    get_block(a,invblock_mm,n,m,i,j);
                    get_block(a,diagblock_mm,n,m,r,j);
                    mat_mult_sub(diagblock_mm,block_mm,invblock_mm,m,m,m);
                    set_block(a,diagblock_mm,n,m,r,j);
                }

                if (is_l!= 0) {
                get_block_ml(a,tmpblock_ml,n,m,l,i);
                get_block_ml(a,tmpblock_ml1,n,m,l,r);
                mat_mult_sub(tmpblock_ml1,block_mm,tmpblock_ml,m,l,m);
                set_block_ml(a,tmpblock_ml1,n,m,l,r);
                }
            }else
            {
                // printlxn(a,n,n,n,n);
               get_block_lm(a, block_ml, n, m, l, i);

            //    printf("block_lm in col %d\n",i);
            //    printlxn(block_ml,m,l,m,n);

               get_block_lm(a, tmpblock_ml, n, m, l, i);
               memset(tmpblock_ml,0,m*l*sizeof(double));
               set_block_lm(a, tmpblock_ml, n, m, l, i);

               get_vec_block(b,vecb_m,n,m,i);// get_vec_block(b,vecb_m,n,m,i);//вычитание из вектора b block_mm*b
            //    cout<<"vecb_m in subtract i= "<<i<<" r="<<r<<endl;
            //     printlxn(vecb_m,m,1,m,m);
               get_vec_block(b,tmpvecb_l,n,m,r);// get_vec_block(b,tmpvecb_m,n,m,r);
               vec_mult_sub_lm(tmpvecb_l,block_ml,vecb_m,l,m);// vec_mult_sub(tmpvecb_m,block_mm,vecb_m,m);
               set_vec_block(b,tmpvecb_l,n,m,r);  // set_vec_block(b,tmpvecb_m,n,m,r);


                for(int j = i + 1; j < k; j++) {
                get_block(a,tmpblock_mm,n,m,i,j);
                get_block_lm(a, tmpblock_ml, n, m, l, j);

                // cout<<"tmpblock_ml in col "<<j<<endl;
                // printlxn(tmpblock_ml,m,l,m,m);

                mat_mult_sub(tmpblock_ml,block_ml,tmpblock_mm,l,m,m);

                
                // get_vec_block(b,vecb_m,n,m,i);//вычитание из вектора b block_mm*b
                // get_vec_block(b,tmpvecb_m,n,m,r);//
                vec_mult_sub_lm(tmpvecb_m,block_ml,vecb_m,l,m);//
                set_block_lm(a, tmpblock_ml, n, m, l, j);
                // set_vec_block(b,tmpvecb_m,n,m,r);//set_vec...
                }

                if (is_l != 0) {
                    get_block_ml(a,tmpblock_ml,n,m,l,i);
                    get_block(a,tmpblock_ll,n,m,k,k);
                    mat_mult_sub(tmpblock_ll,block_ml,tmpblock_ml,l,l,m);
                    set_block(a,tmpblock_ll,n,m,k,k);
                }

            }
            
        }

        // cout<<"LAST PRINT"<<endl;
        // printlxn(a,n,n,n,n);
        // printlxn(b,n,1,n,n);
         
        
    }
       
    
    //начало обратного хода
        

    // cout<<"colsw : "<<endl;
    // for(int i = 0 ; i < k ; i++) cout<<colsw[i]<<" ";


    for(int i = n-1; i >= 0 ; i--)
    {
        if(i == n-1) x[i] = b[i];
        
        else
        {
            x[i] = b[i];
            for(int j = n-1 ; j >i;j--)
            {
                x[i] -= a[i*n + j]*x[j];
            }
        }
    }

    // cout<<"Vector x before swap :"<<endl;
    // printlxn(x,n,1,n,n);

    

    for(int i = 0 ; i < k ; i++)
    {
        

        if(i != colsw[i]){ 
            int t;
            swap_block_vec(x,n,m,i,colsw[i]);
            t = colsw[colsw[i]];
            colsw[colsw[i]] = colsw[i];
            colsw[i] = t; 
        }


    }

    

    return 0;

}






// void pllinit_matrix(double *a,int s, int n , int m , int k, int p)
// {
//     int i,i2,j;
//     pthread_barrier_t b;
//     pthread_barrier_init(&b,0,p);

//     for(i = k*m; i<n; i+=m*p)
//     {
//         int h =(i+m<n ? m:n-i);

//         for(i2 = i; i2<i+h ; i2++)
//         {
//             for(j = 0; j < n ; j++)
//             {
//                 a[i2*n + j] = f(s,n,i2,j);
//             }
//         }
//     }

//     // pthread_barrier_wait(&b);

//     pthread_barrier_destroy(&b);

// }

// void pllinit_vectorb(double *b,double *a, int n , int m , int k, int p)
// {

//     int i,i2,j;
//     pthread_barrier_t bar;
//     pthread_barrier_init(&bar,0,p);

//     for(i = k*m; i<n; i+=m*p)
//     {
//         int h =(i+m<n ? m:n-i);

//         for(i2 = i; i2<i+h ; i2++)
//         {
//             double sumbi = 0;
//     //     for(int k = 0; k <(n-1)/2+1 ; k++)
//     //     {
//     //         sumbi+= a[i*n+2*k];
            
            
//     //     }
        
//     //     b[i] = sumbi;
//             for(j = 0; j <n; j+=2)
//                     {
//                         sumbi+= a[i2*n+j];
//                     }

//                     b[i2] = sumbi;
            
//         }
//     }
    

//     // pthread_barrier_wait(&b);

//     pthread_barrier_destroy(&bar);

// }


// void clear(double *block_mm,double *block_ml,double *block_ll,double *tmpblock_mm,double *tmpblock_ml,double *tmpblock_ml1,double *tmpblock_ll,double *invblock_mm,double *invblock_ll,double *diagblock_mm,double *diaginvblock_mm,double *vecb_m,double *vecb_l,double *tmpvecb_m, double *tmpvecb_l,int *colsw)
// {
//     delete []block_mm ;
//     delete []block_ml ;
//     delete []block_ll ;
//     delete []tmpblock_mm ;
//     delete []tmpblock_ml ;
//     delete []tmpblock_ml1 ;
//     delete []tmpblock_ll ;
//     delete []invblock_mm ;
//     delete []invblock_ll ;
//     delete []diagblock_mm ;
//     delete []diaginvblock_mm ;
//     delete []vecb_m ;
//     delete []vecb_l ;
//     delete []tmpvecb_m ;
//     delete []tmpvecb_l ; 
//     delete []colsw ;
// }



// void* parallelSolve(void* ptr)
// {
//     args *ap = (args*) ptr;

//     double *a = ap->a;
//     double *b = ap->b;
//     double *x = ap->x;
//     int n = ap->n;
//     int m = ap->m;
//     int s = ap->s;
//     int r = ap->r;
//     int thr = ap->thr;
//     int p = ap->p;
//     char* name = ap->name;
//     int *mainblocks = ap->mainblocks;
//     double *minnorms = ap->minnorms;
//     pthread_barrier_t *barrier = ap->barrier;
//     pthread_mutex_t *mutex = ap->mutex;
//     cpu_set_t cpu;
//     static bool isout = false;

//     CPU_ZERO(&cpu);

//     int n_cpus = get_nprocs();
//     int cpu_id = n_cpus -1 -(thr%n_cpus);

//     CPU_SET(cpu_id,&cpu);
//     pthread_t tid = pthread_self();

//     pthread_setaffinity_np(tid,sizeof(cpu),&cpu);

//     //obnulaem matricu chtobi privyazat ee k cpu

//     b=b;

//     // for(int i = thr*m; i<n; i+=p*m)
//     // {
//     //     // int h = (i+m < n ? m : i+m-n);
//     //     // memset(a+i*h,0,h*n*sizeof(double));
//     //     // memset(b+i,0,n*sizeof(double));
//     // }
    

//     if(name)
//     {
//         static int res = 0;
//         res=res;

//         // printf("YA TUT IN THREAD %d",thr);
        

//         if(thr == 0)
//         {
//             res = readarray(a,n,name);
//         }

//         pthread_barrier_wait(barrier);

//         if(res < 0) //tut kazdui thread dolzen znat chto res < 0
//         {
//             if(thr == 0) printf("File %s is bad\n",name);
//             pthread_mutex_lock(mutex);
//             isout = true;
//             pthread_mutex_unlock(mutex);
//         }

    

//     }else
//     {
//         pllinit_matrix(a,s,n,m,thr,p);
        
//     }
//     // init_b(a,s,n,b,thr,p);

//      pthread_barrier_wait(barrier);
//         if(isout)
//         {
//             // printf("ISOUT 1\n");
//             // *ap->r1 = -1;
//             // *ap->r2 = -1;
//             return (void*)-1;
//         }

//     double t = get_cpu_time();


//     pthread_barrier_wait(barrier);

//     if(thr == 0)
//     {
        
//     cout<<"\n MATRIX A :\n";
//     printlxn(a,n,n,n,r);

//     }

    

//     if(normofmatrix(a,n) < EPS64)
//     {
//         printf("Norm of matrix A < 1e-64 \n");
//         // r1 = -1;r2 = -1; // nado perenesti v main
        
//         pthread_mutex_lock(mutex);
//         isout = true;
//         pthread_mutex_unlock(mutex);
//     }


//     pthread_barrier_wait(barrier);
//         if(isout)
//         {
//             // printf("ISOUT 2\n");
//             // *ap->r1 = -1;
//             // *ap->r2 = -1;
//             return (void*)-1;
//         }
//     //print matrix a

//     double eps = 1e-15*normofmatrix(a,n);
    
//     //init for b
//     if(thr == 0) cout<<"\nVector b : \n";

//     pllinit_vectorb(b,a,n,m,thr,p);
//      pthread_barrier_wait(barrier);
    
//     if(thr == 0) printlxn(b,n,1,n,r);


//     //main algorithm

//     int  k, l;

//     k = n/m; l = n - m*k ;
    
//     int is_l = (l == 0) ? 0:1; 

//     double *block_mm = new double[m*m];
//     double *block_ml = new double[m*l];
//     double *block_ll = new double[l*l];
//     double *tmpblock_mm = new double[m*m];
//     double *tmpblock_ml = new double[m*l];
//     double *tmpblock_ml1 = new double[m*l];
//     double *tmpblock_ll = new double[l*l];
//     double *invblock_mm = new double[m*m];
//     double *invblock_ll = new double[l*l];
//     double *diagblock_mm = new double[m*m];
//     double *diaginvblock_mm = new double[m*m];
//     double *vecb_m = new double[m];
//     double *vecb_l = new double[l];
//     double *tmpvecb_m = new double[m];
//     double *tmpvecb_l = new double[l]; 
//     /*static*/ int *colsw = new int[k];

   

//     // int tcol = (k + is_l<p ? 1:(k+is_l)/p);//how many blocks thread have
    

//     for(int c =0; c < k  ;c++) colsw[c]=c;
//     (void)mainblocks;
//     (void)minnorms;
//     (void)mutex;
//     (void)eps;
//     (void)is_l;
//     // (void)mainBlock;

//     pthread_barrier_wait(barrier);

//     for(int i = 0; i < k + is_l; i++)
//     {
//         int mainBlock = -1;
//         double minNorm = 1e64;
//         double locmin = 1e64;
//         int localMainBlock = i;
//         (void)minNorm;
//         mainblocks[thr] = -1;
//         minnorms[thr] = -1;

//         pthread_barrier_wait(barrier);
//         // printf("IN THREAD %d startzone = %d endzone = %d i = %d\n",thr,startzone,endzone,i);
//         if(i != k)
//         {
//             for(int j = i + thr ; j < k ; j+=p)
//             {
//                 get_block(a,block_mm,n,m,i,j);
//                 // printf("Block[%d,%d]\n",i,j);
//                 // printlxn(block_mm,m,m,m,m);

//                 if(inverse(invblock_mm,block_mm,m,eps))
//                 {
//                     if(normofmatrix(invblock_mm,m) < locmin) 
//                     {
                        
//                         locmin = normofmatrix(invblock_mm,m);
//                         localMainBlock = j;
//                         mainblocks[thr] = localMainBlock;
//                         minnorms[thr] = locmin;
                       
//                     }
//                 }
//             }
//         }
//         else
//         {
//             if(thr == 0)
//             {
//                 get_block(a,block_ll,n,m,k,k);

//                 if(!inverse(invblock_ll,block_ll,l,eps))
//                 {
//                     printf("Block [%d,%d] (block[l,l] in our matrix)  has no inverse after the transformations\n",k,k);
//                     pthread_mutex_lock(mutex);
//                     isout = true;
//                     pthread_mutex_unlock(mutex);
//                 }
//                 else
//                 {


                
                
//         //             printlxn(invblock_ll,l,l,l,l);
//                 // minNorm = normofmatrix(invblock_ll,l);
//                 locmin = normofmatrix(invblock_ll,l);
//                 localMainBlock = k;
//                 mainblocks[thr] = localMainBlock;
//                 minnorms[thr] = locmin;
//                 }
//             }   
            
//         }

//             pthread_barrier_wait(barrier);
//             if(isout)
//             {
//                 // printf("ISOUT 3\n");
//                 CLEAR;
//                 // *ap->r1 = -1;
//                 // *ap->r2 = -1;
//                 return (void*)-1;
//             }

//             locmin = 1e64;
//             if(thr == 0)
//             {
//                 for(int q = 0; q < p; q++) // if row have not mainblock
//                 {
//                     int cc = 0;
//                     if(mainblocks[q] == -1)
//                             cc++;
//                     if(cc == p)
//                     {    
//                         printf("Have no main block in row %d\n",i);
//                         // clear(block_mm,block_ml,block_ll,tmpblock_mm,tmpblock_ml,tmpblock_ml1,tmpblock_ll,invblock_mm,invblock_ll,diagblock_mm,diaginvblock_mm,vecb_m,vecb_l,tmpvecb_m, tmpvecb_l,colsw);
//                         pthread_mutex_lock(mutex);
//                         isout = true;
//                         pthread_mutex_unlock(mutex);
                        
//                     }

                    
//                     if(locmin > minnorms[q] && minnorms[q]>=0)
//                         {
//                             locmin = minnorms[q];
//                             mainBlock = mainblocks[q];
//                             minNorm = minnorms[q];
//                             // printf("mainBlock = %d minNorm = %lf thread %d\n",mainBlock,minNorm,thr);
//                         }

//                 }
//             }
//             pthread_barrier_wait(barrier);

//             if(isout)
//             {
//                 CLEAR;
//                 // printf("ISOUT 4\n");
//                 // *ap->r1 = -1;
//                 // *ap->r2 = -1;
//                 return (void*)-1;
//             }
//             pthread_barrier_wait(barrier);

//             // if(thr == 0) 
//             // {    
//             // printf("global mainBlock = %d global minNorm = %lf [i = %d]\n",mainBlock,minNorm,i);
//             // printf("mainBlock :\n");
//             // get_block(a,block_mm,n,m,i,mainBlock);
//             // printlxn(block_mm,m,m,m,r);
//             // printf("invmainBlock :\n");
//             // inverse(invblock_mm,block_mm,m,eps);
//             // printlxn(invblock_mm,m,m,m,r);
//             // printf("block[%d,%d] :\n",1,0);
//             // get_block(a,block_mm,n,m,1,0);
//             // printlxn(block_mm,m,m,m,r);
//             // inverse(invblock_mm,block_mm,m,eps);
//             // printf("block[%d,%d] norm = %lf :\n",1,0,normofmatrix(invblock_mm,m));
//             // printlxn(invblock_mm,m,m,m,r);
//             // }

//             if(thr == 0)
//             {    
//                 if((fabs(minNorm - 1e64) < eps))
//                 {   
//                     if(i!=0)
//                     {
//                         printf("No inverse matrix in row %d after the transformations in thread %d\n",i,thr);
//                         // cout<<"\n MATRIX A :\n";
//                         // printlxn(a,n,n,n,r);
//                         pthread_mutex_lock(mutex);
//                         isout = true;
//                         pthread_mutex_unlock(mutex);
//                     }
//                     else
//                         printf("No inverse matrix in row %d\n",i);

//                     // clear(block_mm,block_ml,block_ll,tmpblock_mm,tmpblock_ml,tmpblock_ml1,tmpblock_ll,invblock_mm,invblock_ll,diagblock_mm,diaginvblock_mm,vecb_m,vecb_l,tmpvecb_m, tmpvecb_l,colsw);

//                         pthread_mutex_lock(mutex);
//                         isout = true;
//                         pthread_mutex_unlock(mutex);
//                 }
//             }

//             pthread_barrier_wait(barrier);
//             if(isout)
//             {
//                 CLEAR;
//                 // printf("ISOUT 5\n");
//                 // *ap->r1 = -1;
//                 // *ap->r2 = -1;
//                 return (void*)-1;
//             }

//             if(mainBlock != i && thr == 0 && i < k)
//             {

//                 // printf("BEFORE SWAP COLUMNS %d %d\n",i,mainBlock);
//                 // printlxn(a,n,n,n,r);
//                 swap_block_columns(a,n,m,i,mainBlock);
//                 // printf("AFTER SWAP COLUMNS %d %d\n",i,mainBlock);
//                 // printlxn(a,n,n,n,r);
//                 swap(colsw[i],colsw[mainBlock]);
//                 cout<<"swapped "<< i<<" "<<mainBlock<<" in row "<<i<<endl;
//             }
//             else if(mainBlock == -1 && thr == 0 && i < k)
//             {
//                 printf("NO mainblock in row %d\n",i);
//                 pthread_mutex_lock(mutex);
//                 isout = true;
//                 pthread_mutex_unlock(mutex);
//             }

            
//             pthread_barrier_wait(barrier);
//             if(isout)
//             {
//                 // printf("CHECK0 2601str THREAD %d\n",thr);
//                 CLEAR;
//                 // printf("ISOUT 6\n");
//                 // *ap->r1 = -1;
//                 // *ap->r2 = -1;
//                 return (void*)-1;
//             }
//             pthread_barrier_wait(barrier);//important

//             // printf("BEFORE MULT thread %d\n",thr);
//             // printlxn(a,n,n,n,r);

//             //start multiplication
//         if(i<k)
//         {
//             get_block(a,diagblock_mm,n,m,i,i);
//             // printf("CHECK1 THREAD %d\n",thr);

//             if(!(inverse(diaginvblock_mm,diagblock_mm,m,eps)))
//             {
//                 printf("no blocks in row has inverse block i = %d thread %d\n",i,thr);
//                 // printlxn(a,n,n,n,r);
//                 // CLEAR;
//                 pthread_mutex_lock(mutex);
//                 isout = true;
//                 pthread_mutex_unlock(mutex);
//             }

//             pthread_barrier_wait(barrier);
//             if(isout)
//             {
//                 CLEAR;
//                 // printf("ISOUT 7\n");
//                 // *ap->r1 = -1;
//                 // *ap->r2 = -1;
//                 return (void*)-1;
//             }

//             // printf("CHECK2 THREAD %d\n",thr);

            
       
//             if(thr == 0)
//             {
//                 get_vec_block(b,vecb_m,n,m,i);
//                 mat_x_vector(tmpvecb_m,diaginvblock_mm,vecb_m,m);// double *resvec = mat_x_vector(diaginvblock_mm,vecb_m,m);
//                 // cout<<"tmpvecb_m : "<<endl;
//                 // printlxn(tmpvecb_m,m,1,m,m);    
//                 set_vec_block(b,tmpvecb_m,n,m,i);
//             }

            

            

//             for(int j = thr + i ; j < k ; j+=p) //mb try j = i
//             {
//                 // pthread_mutex_lock(mutex);
//                 get_block(a,block_mm,n,m,i,j);
                
//                 multiplication(tmpblock_mm,diaginvblock_mm,block_mm,m,m,m);// matmult(tmpblock_mm,diaginvblock_mm,block_mm,m,m,m);// double *resmult = matmult(diaginvblock_mm,block_mm,m,m,m)

//                 set_block(a,tmpblock_mm,n,m,i,j);
                
//                 if (!block_mm || !vecb_m || !invblock_mm) 
//                 {
//                     fprintf(stderr, "Error: temporary buffers not initialized!\n");
//                     // CLEAR;
//                     // return (void*)-1;
//                     pthread_mutex_lock(mutex);
//                     isout = true;
//                     pthread_mutex_unlock(mutex);
//                 }
//                 // pthread_mutex_unlock(mutex);
//             }

            
//             if(is_l != 0 && thr == p-1)
//             {
//                 get_block_ml(a,block_ml,n,m,l,i);
//                 multiplication(tmpblock_ml,diaginvblock_mm,block_ml,m,m,l);// matmult(tmpblock_ml,diaginvblock_mm,block_ml,m,m,l);
//                 set_block_ml(a,tmpblock_ml,n,m,l,i);
//             }

//             pthread_barrier_wait(barrier);

//             if(isout)
//             {
//                 CLEAR;
//                 // printf("ISOUT 8\n");
//                 // *ap->r1 = -1;
//                 // *ap->r2 = -1;
//                 return (void*)-1;
//             }
            

//             // if(thr == 0 && i == 0) 
//             //     {
//             //         printf("Matrix A after mult \n");
//             //         printlxn(a,n,n,n,r);
//             //     }
        
//         }
//         else //if(i == k && thr == 0) uncomm later when i do subtract 
//         {   
//             // printlxn(a,n,n,n,n);
//             // printlxn(b,n,1,n,n);
//             if(thr == p-1)
//             {           
//             // pthread_barrier_wait(barrier);
//             get_block(a,block_ll,n,m,i,i);
//             get_vec_block(b,vecb_l,n,m,i);

//             // printf("Block ll :\n");
//             // printlxn(block_ll,l,l,l,n);

//             // cout<<"vecb_l:"<<endl;
//             // printlxn(vecb_l,l,1,l,n);

//             if(!(inverse(invblock_ll,block_ll,l,eps)))
//                 {
//                     printf("ll block has no inverse in thread %d\n",thr);
//                     // CLEAR; // нужно сделать не выход ретурн -1 а завести флаг по которому потом выйдут одновременно все потоки а то будет DL
//                     pthread_mutex_lock(mutex);
//                     isout = true;
//                     pthread_mutex_unlock(mutex);
//                 }
//                 // pthread_barrier_wait(barrier);
            
//             // printlxn(invblock_ll,l,l,l,n);

//             multiplication(tmpblock_ll,invblock_ll,block_ll,l,l,l);// matmult(tmpblock_ll,invblock_ll,block_ll,l,l,l);


//             // for(int ii = 0 ; ii < l; ii++)
//             // {
//             //     for(int jj = 0 ; jj < l; jj++)
//             //     {
//             //         if(ii != jj) block_ll[ii*l+jj] = 0;
//             //         else block_ll[ii*l+jj] = 1;
//             //     }
//             // }

//             mat_x_vector(tmpvecb_l,invblock_ll,vecb_l,l);
            
//             // printlxn(tmpvecb_l,l,1,l,n);

//             set_block(a,tmpblock_ll,n,m,i,i);
//             set_vec_block(b,tmpvecb_l,n,m,i);
//             }
//             // cout<<"WE ARE IN i = k"<<endl;
//         }

//         pthread_barrier_wait(barrier);

//             if(isout)
//             {
//                 CLEAR;
//                 // printf("ISOUT 9\n");
//                 // *ap->r1 = -1;
//                 // *ap->r2 = -1;
//                 return (void*)-1;
//             }

//         // if(thr == 0 && i == 0) 
//         //         {
//         //             printf("Matrix A before subtr \n");
//         //             printlxn(a,n,n,n,r);
//         //         }
//         pthread_barrier_wait(barrier);

//         for(int rr = thr + i + 1 ; rr < k + is_l ; rr+=p)
//         {
//             if(rr < k)
//             {  
//                 get_block(a,block_mm,n,m,rr,i);
                

//                 // pthread_barrier_wait(barrier);

//                 // pthread_mutex_lock(mutex);
//                 get_block(a,tmpblock_mm,n,m,rr,i);
//                 // printf("tmpblock_mm before MEMSET:\n");
//                 // printlxn(tmpblock_mm,m,m,m,r);
//                 memset(tmpblock_mm,0, m*m*sizeof(double));
//                 set_block(a,tmpblock_mm,n,m,rr,i);
//                 // printf("A after MEMSET:\n");
//                 // printlxn(a,n,n,n,r);
//                 // pthread_mutex_unlock(mutex);

//                 // pthread_barrier_wait(barrier);

//                 // pthread_mutex_lock(mutex); // mb delete need test
                

//                 // if(thr == 0)
//                 // {
                    

//                     // not in i for
//                     get_vec_block(b,vecb_m,n,m,i);//вычитание из вектора b block_mm*b
//                     get_vec_block(b,tmpvecb_m,n,m,rr);
//                     vec_mult_sub(tmpvecb_m,block_mm,vecb_m,m);
//                     set_vec_block(b,tmpvecb_m,n,m,rr);
//                 // }

//                 // cout<<"tmpvecb_m in subtract i= "<<i<<" r="<<r<<endl;
//                 // printlxn(tmpvecb_m,m,1,m,m);
//                 // printlxn(b,n,1,n,n);

//                 for (int j = i + 1; j < k; j++) { //(int j = i + 1; j < k; j++)
                    
//                     get_block(a,invblock_mm,n,m,i,j);
//                     get_block(a,diagblock_mm,n,m,rr,j);
//                     mat_mult_sub(diagblock_mm,block_mm,invblock_mm,m,m,m);
//                     set_block(a,diagblock_mm,n,m,rr,j);
                    
//                 }
//                 // pthread_mutex_unlock(mutex);

//                 if (is_l!= 0) {
//                 get_block_ml(a,tmpblock_ml,n,m,l,i);
//                 get_block_ml(a,tmpblock_ml1,n,m,l,rr);
//                 mat_mult_sub(tmpblock_ml1,block_mm,tmpblock_ml,m,l,m);
//                 set_block_ml(a,tmpblock_ml1,n,m,l,rr);
//                 }
//             }
//             else //if(rr == k && thr == 0)
//             {
//                 // printlxn(a,n,n,n,n);
//                 // if(thr == 0)
//                 // {
//                 get_block_lm(a, block_ml, n, m, l, i);

//             //    printf("block_lm in col %d\n",i);
//             //    printlxn(block_ml,m,l,m,n);

//                get_block_lm(a, tmpblock_ml, n, m, l, i);
//                memset(tmpblock_ml,0,m*l*sizeof(double));
//                set_block_lm(a, tmpblock_ml, n, m, l, i);

//             //    pthread_barrier_wait(barrier);

//                get_vec_block(b,vecb_m,n,m,i);// get_vec_block(b,vecb_m,n,m,i);//вычитание из вектора b block_mm*b
//             //    cout<<"vecb_m in subtract i= "<<i<<" r="<<r<<endl;
//             //     printlxn(vecb_m,m,1,m,m);
//                get_vec_block(b,tmpvecb_l,n,m,rr);// get_vec_block(b,tmpvecb_m,n,m,r);
//                vec_mult_sub_lm(tmpvecb_l,block_ml,vecb_m,l,m);// vec_mult_sub(tmpvecb_m,block_mm,vecb_m,m);
//                set_vec_block(b,tmpvecb_l,n,m,rr);  // set_vec_block(b,tmpvecb_m,n,m,r);


//                 for(int j = i + 1; j < k; j++) {
//                 get_block(a,tmpblock_mm,n,m,i,j);
//                 get_block_lm(a, tmpblock_ml, n, m, l, j);

//                 // cout<<"tmpblock_ml in col "<<j<<endl;
//                 // printlxn(tmpblock_ml,m,l,m,m);

//                 mat_mult_sub(tmpblock_ml,block_ml,tmpblock_mm,l,m,m);

                
// //                 get_vec_block(b,vecb_m,n,m,i);//вычитание из вектора b block_mm*b
//                 get_vec_block(b,tmpvecb_m,n,m,rr);//
//                 vec_mult_sub_lm(tmpvecb_m,block_ml,vecb_m,l,m);//
//                 set_block_lm(a, tmpblock_ml, n, m, l, j);
//                 // }
//                 // set_vec_block(b,tmpvecb_m,n,m,r);//set_vec...
//                 }

//                 if (is_l != 0) {
//                     get_block_ml(a,tmpblock_ml,n,m,l,i);
//                     get_block(a,tmpblock_ll,n,m,k,k);
//                     mat_mult_sub(tmpblock_ll,block_ml,tmpblock_ml,l,l,m);
//                     set_block(a,tmpblock_ll,n,m,k,k);
//                 }

//             }
//             // pthread_barrier_wait(barrier);
//             // if(thr == 0)
//             // {
//             //     printf("Matrix A after subtr: \n");
//             //     printlxn(a,n,n,n,r);
//             // }
//         }
//     }

//     pthread_barrier_wait(barrier);

//             if(isout)
//             {
//                 CLEAR;
//                 // printf("ISOUT 10\n");
//                 // *ap->r1 = -1;
//                 // *ap->r2 = -1;
//                 return (void*)-1;
//             }

//     // if(thr == 0)
//     // {
//     //     printf("Matrix A AFTER ALL OPERATIONS: \n");
//     //     printlxn(a,n,n,n,r);
//     //     printf("Vector b AFTER ALL OPERATIONS: \n");
//     //     printlxn(b,n,1,n,r);
//     // }

    
//     //start reverse alg
        
//     if(thr == 0)
//     {
//             for(int i = n-1; i >= 0 ; i--)
//         {
//             if(i == n-1) x[i] = b[i];
            
//             else
//             {
//                 x[i] = b[i];
//                 for(int j = n-1 ; j >i;j--)
//                 {
//                     x[i] -= a[i*n + j]*x[j];
//                 }
//             }
//         }

//         // cout<<"Vector x before swap :"<<endl;
//         // printlxn(x,n,1,n,n);

        

//         for(int i = 0 ; i < k ; i++)
//         {
            

//             if(i != colsw[i]){ 
//                 int t;
//                 swap_block_vec(x,n,m,i,colsw[i]);
//                 t = colsw[colsw[i]];
//                 colsw[colsw[i]] = colsw[i];
//                 colsw[i] = t; 
//             }


//         }
//     }
//     pthread_barrier_wait(barrier);
//     // tut osvobozdaem memory allocated in thread
//     CLEAR;

//     if(name)
//     {

//         // printf("YA TUT IN THREAD %d",thr);
       
//         if(thr == 0)
//         {
//             readarray(a,n,name);
//         }

//     }else
//     {
//         pllinit_matrix(a,s,n,m,thr,p);
//         // pthread_barrier_wait(barrier);
//     }

//     pthread_barrier_wait(barrier);
//         // if(isout)
//         // {
//         //     return (void*)-1;
//         // }
//     // tut osvobozdaem memory allocated in thread

//     // if(thr == 0)
//     // {
        
//     // cout<<"\n MATRIX A REINIT :\n";
//     // printlxn(a,n,n,n,r);

//     // }
   

//     // if(thr == 0) cout<<"\nVector b REINIT : \n";

//     pllinit_vectorb(b,a,n,m,thr,p);
    
//     // if(thr == 0) printlxn(b,n,1,n,r);


//     t = get_cpu_time() - t;

//     ap->time =t;

//     printf("CPU Time thread %d = %.2lf\n",thr,ap->time);

//     pthread_barrier_wait(barrier);

//     if(thr == 0)
//     {
//         cout<<"\nSolution vector x : ";
//         printlxn(x,n,1,n,r);
//     }

//     return (void*)1;
// }


int l2g(/*int n,*/ int m, int k, int p, int i_loc)
{
    int i_loc_m = i_loc/m;
    int i_glob_m = i_loc_m*p + k;
    return i_glob_m*m +  i_loc % m;
}

int g2l(/*int n,*/ int m, /*int k,*/ int p, int i_glob)
{
    int i_glob_m = i_glob/m;
    int i_loc_m = i_glob_m%p;
    return i_loc_m*m + i_glob%m;
}

int l2g_block(/*int n, int m,*/ int k, int p, int i_loc_m)
{
    // int i_loc_m = i_loc/m;
    return  i_loc_m*p + k;

}

int g2l_block(/*int n, int m, int k,*/ int p, int i_glob_m)
{
    // int i_glob_m = i_glob/m;
    return i_glob_m%p;
    
}

int get_k(/*int n,*/ int m, int p, int i_glob)
{
    int i_glob_m = i_glob/m;
    return i_glob_m%p;
}

int get_block_rows(int n, int m, int p, int k)
{
    int b = (n + m - 1)/m;
    return (b%p <= k ? b/p : b/p + 1);
}

int get_max_block_rows(int n, int m, int p)
{
    int b = (n + m - 1)/m;
    return (b + p - 1)/p;
}
int get_rows(int n, int m, int p, int k)
{
    int b = (n + m - 1)/m;
    int b_last = (b - 1)%p;
    int b_loc = (b%p <= k ? b/p : b/p + 1); // num of local block rows

    if(k != b_last) return b_loc * m;

    int kk = n/m;
    int l = n - kk*m;

    if(l == 0) return b_loc * m;
    if(b_loc == 0) return 0;

    return (b_loc - 1)*m + l;
}

void init_matrix(double *a, int n, int m, int p, int k, double (*f)(int s,int n,int i,int j),int s)
{
    int i_loc, i_glob, j_loc, j_glob, rows;
    rows = get_rows(n,m,p,k);
    for(i_loc = 0; i_loc < rows; i_loc++)
    {
        i_glob = l2g(m,k,p,i_loc);
        for(j_loc = 0; j_loc < n; j_loc++)
        {
            j_glob = j_loc;
            a[i_loc*n + j_loc] = f(s,n,i_glob,j_glob);
        }
    }
}

void init_vector_b(double *a,double *b,int n, int m, int p, int k)
{
    int i_loc,j_loc,rows;
    rows = get_rows(n,m,p,k);
    memset(b,0,rows*sizeof(double));

    for(i_loc = 0; i_loc < rows; i_loc++)
    {
        double sum = 0;
        for(j_loc = 0; j_loc < n; j_loc++)
        {
            sum += a[i_loc * n + j_loc] * ((j_loc+1)%2);
            // printf("a[i_loc * n + j_loc] = %10.3e\n",a[i_loc * n + j_loc]);
        }
        b[i_loc] = sum;
    }
    // printf("Proc %d ",k);
    // for(i_loc = 0; i_loc < rows; i_loc++)
    //     printf("%10.3e ",b[i_loc]);
    // printf("\n");
}

int read_matrix(double *a, int n, int m, int p, int k, const char *name,double *buf, MPI_Comm com)
{
    int main_k = 0;
    FILE *fp = nullptr; int err = 0;
    // if (k == main_k) printf("Name in read_matrix %s\n",name);

    if(k == main_k)
    {
        fp = fopen(name,"r");
        if(!fp) 
        {
            err = 1;
            printf("have error: bad file \n");
        }
    }
    MPI_Bcast(&err,1,MPI_INT,main_k,com);

    if(err) return err;

    memset(buf,0,n*m*sizeof(double));

    int b, max_b = (n + m -1)/m;
    int printed = 0;

    for(b = 0; b < max_b; b++)
    {
        int owner = b%p; //whos row now
        int rows = (b*m + m <=n ? m: n - b*m);// rows in block row
        int b_loc = b/p; // local block num
        

        if(k == main_k)
        {
            err += read_array(fp,buf,printed,n,n*rows);
            // if(printed > n*n)
            // {
            //     printf("error size of file\n");
            //     return -1;
            // }
            
            // for(int i = 0; i < rows; i++)
            // {
            //     for(int j = 0; j < n; j++)
            //         {
            //             printf("%10.3e ",buf[i*n+j]);
            //         }
            //         printf("\n");
            // }
            if(owner == main_k)
                memcpy(a + b_loc*m*n,buf,n*rows*sizeof(double));
            else
                MPI_Send(buf,n*rows,MPI_DOUBLE,owner,0,com);
        }
        else if (owner == k)
        {
            MPI_Status st;
            MPI_Recv(a + b_loc*n*m,n*rows,MPI_DOUBLE,main_k,0,com,&st);
        }

        
    }

    if(k == main_k)
    {
        fclose(fp);
        fp = nullptr;
        if(printed != n*n) err = 1;
    }
    MPI_Bcast(&err,1,MPI_INT,main_k,com);
    if(err) 
    {
        if(k == main_k) printf("error\n");
        return err;
    }

    return 0;
}
void print_matrix(double *a, int n, int m, int p, int k, double *buf/*nxm size*/, int max_print, MPI_Comm com)
{
    int main_k = 0;
    int b, max_b = (n + m -1)/m;
    int printed_rows = 0;

    for(b = 0; b < max_b; b++)
    {
        int owner = b%p, b_loc = b/p;
        int rows = min(m,n - b*m);

        if(k == main_k)
        {
            if(owner == main_k)
                printed_rows += print_array(a + b_loc*n*m,n,rows,printed_rows,max_print);
            else
            {
                MPI_Status st;
                MPI_Recv(buf, n*rows, MPI_DOUBLE, owner,0,com,&st);
                printed_rows += print_array(buf,n,rows,printed_rows,max_print);
            }
        }
            else if(k == owner)
                MPI_Send(a + b_loc*m*n, n*rows,MPI_DOUBLE,main_k,0,com);
        
    }
}

void print_vector(double *b,int n,int m, int p, int k,double *buf,int max_print,MPI_Comm com)
{
    int main_k = 0;
    int i_loc/*i_glob,i_glob_m,rows*/;
    int bl, max_b = (n + m -1)/m;
    int r = min(n,max_print);
    int printed = 0;
    

    for(bl = 0; bl < max_b; bl++)
    {
        int owner = bl%p, b_loc = bl/p;
        b_loc=b_loc;
        int rows = min(m,n - bl*m);

        if(k == main_k)
        {
            if(owner == main_k)
                {
                    for(i_loc = 0; (printed < r)&&(i_loc < rows);printed++,i_loc++)
                        printf("%10.3e ",b[b_loc*m + i_loc]);
                }
            else
            {
                MPI_Status st;
                MPI_Recv(buf,rows, MPI_DOUBLE, owner,0,com,&st);
                for(i_loc = 0; (printed < r)&&(i_loc < rows);printed++,i_loc++)
                        printf("%10.3e ",buf[i_loc]);
            }
        }
        else if(k == owner)
        {
            MPI_Send(b + b_loc*m,rows,MPI_DOUBLE,main_k,0,com);
        }
        
    }
}

int print_array(double *a, int n, int m, int printed_rows,int max_print)
{
    if(printed_rows >= max_print) return 0;
    int p_n = (n > max_print ? max_print : n);
    int p_m = (printed_rows + m <= max_print ? m:max_print - printed_rows);

    for(int i = 0; i < p_m; i++)
    {
        for(int j = 0; j < p_n; j++)
        {
            printf(" %10.3e", a[i*n + j]);
        }
        printf("\n");
    }
    return p_m;
}

void matrix_mult_vector(double *a, double *bvec,double *c, int n, int m, int k, int p,MPI_Comm com)
{
    int b_rows = get_block_rows(n,m,p,k);
    int rows = get_rows(n,m,p,k);
    int max_b_rows = get_max_block_rows(n,m,p);
    // int b, max_b = (n + m - 1)/m;
    int dst = (k-1+p)%p;
    int src = (k+1+p)%p;

    memset(c,0,rows*sizeof(double));

    for(int s = 0; s < p; s++)//resend loop
    {
        int sk = (k+s)%p; // whos data now in vector b
        int sk_rows = get_rows(n,m,p,sk);
        int sk_b_rows = get_block_rows(n,m,p,sk);
        sk_rows=sk_rows;

        for(int i = 0; i < b_rows; i++)//loop po svoim strokam
        {
            int i_glob_m = l2g_block(k,p,i);
            int h = (i_glob_m*m + m < n ? m:n-i_glob_m*m);

            for(int sk_i = 0; sk_i < sk_b_rows; sk_i++)
            {
                int i_s_glob = l2g_block(sk,p,sk_i);
                int w = (i_s_glob*m + m < n ? m:n-i_s_glob*m);

                for(int ii = 0; ii < h; ii++)
                {
                    double sum = 0;
                    for(int jj = 0; jj < w; jj++)
                    {
                        sum += a[(i*m+ii)*n + i_s_glob*m + jj]*bvec[sk_i*m + jj];
                    }
                    c[i*m + ii] += sum;
                }
            }
        }
        MPI_Status st;
        MPI_Sendrecv_replace(bvec,max_b_rows*m,MPI_DOUBLE,dst,0,src,0,com,&st);

    }
}

void get_block_row(double *a,double *buf,int n, int m, int p,int kk,int i_loc_m)
{
    memset(buf,0,m*n*sizeof(double));

    double *start = a + i_loc_m*m*n;
    int last = get_block_rows(n,m,p,kk)-1;
    int lastnum = l2g_block(kk,p,last);
    int b = n/m;
    int p_m = (lastnum == l2g_block(kk,p,i_loc_m) ? n-b*m:m);

    if(p_m != m)
        printf("in get_block_row p_m = %d proc %d owner = %d i_loc_m = %d\n",p_m,kk,l2g_block(kk,p,i_loc_m)%p,i_loc_m);
    

    for(int i = 0; i < p_m; i++)
        for(int j = 0; j < n; j++)
        {
            buf[i*n+j] = start[i*n+j];
        }

}

int MPI_Solve(double *a, double *b, double *x,int n,int m,int p,int kk,
    double *buf,double *vecbuf,
    double *block_mm, double *block_ml, double *block_ll,
    double *invblock_mm, double *diaginvblock_mm, 
    double *invblock_ll,double *diagblock_mm,
    int *colsw,double *vecb_m,double *vecb_l,
    double *tmpblock_mm,double *tmpblock_ml,double *tmpblock_ml1,double *tmpblock_ll,
    double *tmpvecb_m,double *tmpvecb_l,
    MPI_Comm com)
{
    int main_kk = 0;
    main_kk=main_kk;
    b=b;
    x=x;
    buf=buf;
    vecbuf=vecbuf;
    block_mm=block_mm;
    block_ml=block_ml;
    block_ll=block_ll;
    invblock_mm=invblock_mm;
    diaginvblock_mm=diaginvblock_mm;
    invblock_ll=invblock_ll;
    diagblock_mm=diagblock_mm;
    colsw=colsw;
    vecb_m=vecb_m;
    vecb_l=vecb_l;
    tmpblock_mm=tmpblock_mm;
    tmpblock_ml=tmpblock_ml;
    tmpblock_ml1=tmpblock_ml1;
    tmpblock_ll=tmpblock_ll;
    tmpvecb_m=tmpvecb_m;
    tmpvecb_l=tmpvecb_l;
    


    double N{};
    mpi_norm_matrix(a,N,n,m,p,kk,com);
    if(N < EPS64)
    {
        if(kk == main_kk) printf("Norm of matrix A < 1e-64 \n");
        return -2;
    }
    
    // if(kk == main_kk) printf("Norm of matrix A: %10.3e\n",N);

    int k = n/m, l = n - k*m;
    l=l;
    int is_l = (l==0 ? 0:1);
    memset(colsw,0,k*sizeof(int));
    int b_loc_rows = get_block_rows(n,m,p,kk);
    int owner;
    DI main_block{};
    (void)b_loc_rows;
    main_block=main_block;

    for(int i_glob_m = 0; i_glob_m < k+is_l; i_glob_m++)
    {
        owner = i_glob_m%p;
        int i_loc_m = g2l_block(p,i_glob_m);
        int p_m = (i_glob_m == k ? l:m);
        
        get_block_row(a,buf,n,m,p,owner,i_loc_m);
        if(kk == owner)
        {
            printf("Owner %d printing buf:\n",owner);
            for(int i = 0; i < p_m; i++)
            {
                for(int j = 0; j < n; j++)
                {
                printf(" %10.3e",buf[i*n+j]);
                }
                printf("\n");
            }
        }
        MPI_Bcast(buf,n*p_m,MPI_DOUBLE,owner,MPI_COMM_WORLD);//send block row to all
        // if(kk != owner)
        // {
        //     printf("Printing buf after send proc %d owner = %d:\n",kk,owner);
        //     for(int i = 0; i < p_m; i++)
        //     {
        //         for(int j = 0; j < n; j++)
        //         {
        //         printf(" %10.3e",buf[i*n+j]);
        //         }
        //         printf("\n");
        //     }
        // }
        
        // else
        // {
        //     MPI_Status st;
        //     MPI_Recv(buf,n*p_m,MPI_DOUBLE,owner,0,MPI_COMM_WORLD,&st);
        // }
        main_block.norm = 1e64;
        main_block.num = -1;
    }
    

    return 0;
}