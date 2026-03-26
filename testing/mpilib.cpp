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







void mpi_residuals(double &r1,double &r2,double *a,double *b,double *x,double *realx,int n,int m,int p, int k, double *Ax,MPI_Comm com)
{
    int lrows = get_rows(n,m,p,k);
    memset(Ax,0,lrows*sizeof(double));
    matrix_mult_vector(a,x,Ax,n,m,k,p,com);// mat_x_vector(Ax,a,x,n);// double *Ax = mat_x_vector(a,x,n);

    double ln1{},ld1{},ln2{},ld2{},n1{},d1{},n2{},d2{};

    

    for(int i =0; i < lrows;i++)
    {
        ln1 += fabs(Ax[i]-b[i]);
        ld1 += fabs(b[i]);

        ln2 += fabs(x[i]-realx[i]);
        ld2 += fabs(realx[i]);
    }

    MPI_Allreduce(&ln1,&n1,1,MPI_DOUBLE,MPI_SUM,com);
    MPI_Allreduce(&ln2,&n2,1,MPI_DOUBLE,MPI_SUM,com);
    MPI_Allreduce(&ld1,&d1,1,MPI_DOUBLE,MPI_SUM,com);
    MPI_Allreduce(&ld2,&d2,1,MPI_DOUBLE,MPI_SUM,com);

    r1 = (d1 > 0) ? n1/d1:0;
    r2 = (d2 > 0) ? n2/d2:0;
 
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

 inline void get_block(double *a, double *b, int n, int m, int i_loc_m, int j_loc_m,int i_glob_m)
{
    int i1=0, j1=0, k, l, r, h;
    if(m == 0) {
        printf("m == 0 in i = %d , j = %d\n",i_glob_m,j_loc_m);
        return;
    }
    k = n/m; l = n - m*k ;
    if(i_glob_m < k) r = m;
    else r = l;

    if(j_loc_m < k) h = m;
    else h = l;


    double *bl = a + i_loc_m*n*m + j_loc_m*m; //start of block double *bl = a + i*n*m + j*m;

    for (i1 =0 ; i1 < r ; i1++)
    {
        for (j1=0 ; j1 < h ; j1++)
        {
            b[i1*h + j1] = bl[i1*n + j1];
            

        }
    }
    

}



inline void set_block(double *a, double *b, int n, int m, int i_loc_m, int j_loc_m,int i_glob_m)
{
    int i1=0, j1=0, k, l, r, h;
    k = n/m; l = n - m*k ;
    if(i_glob_m < k) r = m;
    else r = l;

    if(j_loc_m < k) h = m;
    else h = l;


    double *bl = a + i_loc_m*n*m + j_loc_m*m; //start of block

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

double* inverse(double *result,double* A,double *E, int size,double eps)
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

    // double* a = new double[size*size];
    memcpy(result, A, size*size * sizeof(double));


    // double* E = new double[size*size];
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
        double maxEl = fabs(result[i*size+i]);
        int pivot = i;
        for(int j = i+1; j<size;j++)
        {
            if(fabs(result[i*size +j]) > maxEl)
            {
                maxEl = fabs(result[i*size +j]);
                pivot = j;
            }
        }

        
        //swap columns
        if(pivot!=i)
        {
            for(int k = 0; k<size ; k++)
            {
                swap(result[k*size+i],result[k*size+pivot]);
                swap(E[k*size+i],E[k*size+pivot]);
            }
            swap(colsw[i],colsw[pivot]);
        }

        if(fabs(result[i*size+i]) < eps)
        {
            // printf("matrix has no inverse\n");
            // delete []E;
            delete []colsw;
            // delete []a;
            return nullptr;
        }

        //devide row i 
        // cout<<"A matr before devideing"<<endl;
        // printlxn(a,size,size,size,size);

        double mainEl = result[i*size+i];
        for(int k = 0 ; k<size; k++)
        {
            result[i*size+k] /= mainEl;
            E[i*size+k] /= mainEl;
        }
        // cout<<"A matr after devideing"<<endl;
        // printlxn(a,size,size,size,size);


        for(int k = 0; k< size;k++)
        {
            if(k!=i)
            {
                double factor = result[k*size+i];
                for(int j = 0; j <size;j++)
                {
                    result[k*size+j] -= factor * result[i*size+j];
                    E[k*size+j] -= factor * E[i*size+j];
                }
            }
        }

    }

// cout<<"Last a presentation before swap"<<endl;
// printlxn(a,size,size,size,size);


    // swap back columns
    

    // for (int i = 0; i < size*size; ++i) result[i] = 0.0;

    
    
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
            result[ colsw[i]*size + colsw[j] ] = E[i*size + j];
// cout<<"Last a presentation after swap"<<endl;

// printlxn(a,size,size,size,size);

    // delete []a;
    delete [] colsw;
    // delete []E;
    return result;
    
    
}

int inverse1(int m, double *a, double *inv_a, double *tmp_a) {
    
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            inv_a[i * m + j] = (i == j) ? 1.0 : 0.0;
            tmp_a[i * m + j] = a[i * m + j];
        }
    }

   
    for (int i = 0; i < m; ++i) {
        
        int sel = i;
        double max_val = std::abs(tmp_a[i * m + i]);
        for (int k = i + 1; k < m; ++k) {
            double val = std::abs(tmp_a[k * m + i]);
            if (val > max_val) {
                max_val = val;
                sel = k;
            }
        }

        if (max_val < 1e-100) return -1;

     
        if (sel != i) {
            for (int j = 0; j < m; ++j) {
                std::swap(tmp_a[i * m + j], tmp_a[sel * m + j]);
                std::swap(inv_a[i * m + j], inv_a[sel * m + j]);
            }
        }

       
        double diag = 1.0 / tmp_a[i * m + i];
        for (int j = i; j < m; ++j) tmp_a[i * m + j] *= diag;
        for (int j = 0; j < m; ++j) inv_a[i * m + j] *= diag;

     
        for (int k = 0; k < m; ++k) {
            if (k != i) {
                double factor = tmp_a[k * m + i];
  

                for (int j = 0; j < m; ++j) {
                    if (j >= i) tmp_a[k * m + j] -= factor * tmp_a[i * m + j];
                    inv_a[k * m + j] -= factor * inv_a[i * m + j];
                }
            }
        }
    }
    return 0;
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


void swap_block_columns(double *a, int n,int m,int rows, int i, int j)
{
    for(int p =0 ; p < m ; p++)
                {
                    for(int c = 0; c<rows ; c++)
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

inline void get_vec_block(double *b,double *block,int n, int m, int i_loc_m,int i_glob_m)
{
    int i1=0, k, l, r;
    k = n/m; l = n - m*k ;
    if(i_glob_m < k) r = m;
    else r = l;

    for (i1 =0 ; i1 < r ; i1++)
    {
            block[i1] = b[i_loc_m*m+i1];
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

inline void set_vec_block(double *b,double *block,int n, int m, int i_loc_m,int i_glob_m)
{
    int i1=0, k, l, r;
    k = n/m; l = n - m*k ;
    if(i_glob_m < k) r = m;
    else r = l;

    


    //start of block 

    for (i1 =0 ; i1 < r ; i1++)
    {
            b[i_loc_m*m+i1]=block[i1] ;
    }
    

}

inline void get_block_ml(double *a, double *b, int n, int m,int l, int i)
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

inline void set_block_ml(double *a, double *b, int n, int m,int l, int i)
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

inline void get_block_lm(double *a, double *b,int locrows, int n, int m,int l, int j)
{
int i1=0, j1=0;
   

    double *bl = a + (locrows-l)*n + j*m; //start of block double *bl = a + i*n*m + j*m;

    for (i1 =0 ; i1 < l ; i1++)
    {
        for (j1=0 ; j1 < m ; j1++)
        {
            b[i1*m + j1] = bl[i1*n + j1];
            

        }
    }
    
}

inline void set_block_lm(double *a, double *b,int locrows, int n, int m,int l, int j)
{
int i1=0, j1=0;
   

    double *bl = a + (locrows-l)*n + j*m; //start of block double *bl = a + i*n*m + j*m;

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







inline void multiplication(double* Result, double* Block_A, double* Block_B, const int row_A, const int col_row,const int col_B)
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
    return i_glob_m/p;
    
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

void get_block_row(double *a,double *buf,int n, int m, int p_m, int p,int kk,int i_loc_m,int i_glob_m)
{
    // int i_glob_m = l2g_block(kk,p,i_loc_m);
    int owner = i_glob_m%p;
    memset(buf,0,p_m*n*sizeof(double));
    if(owner != kk) return;

    double *start = a + i_loc_m*m*n;
    // int last = get_block_rows(n,p_m,p,kk)-1;
    // int lastnum = l2g_block(kk,p,last);
    // int b = n/m;
    // int p_m = (lastnum == l2g_block(kk,p,i_loc_m) ? n-b*m:m);

    // if(p_m != m)
    //     printf("in get_block_row p_m = %d proc %d owner = %d i_loc_m = %d\n",p_m,kk,l2g_block(kk,p,i_loc_m)%p,i_loc_m);
    

    for(int i = 0; i < p_m; i++)
        for(int j = 0; j < n; j++)
        {
            buf[i*n+j] = start[i*n+j];
        }

}

void set_block_row(double *a,double *buf,int n, int m, int p_m, int p,int kk,int i_loc_m,int i_glob_m)
{
    // int i_glob_m = l2g_block(kk,p,i_loc_m);
    int owner = i_glob_m%p;
    
    if(owner != kk) return;

    double *start = a + i_loc_m*m*n;
    // int last = get_block_rows(n,p_m,p,kk)-1;
    // int lastnum = l2g_block(kk,p,last);
    // int b = n/m;
    // int p_m = (lastnum == l2g_block(kk,p,i_loc_m) ? n-b*m:m);

    // if(p_m != m)
    //     printf("in get_block_row p_m = %d proc %d owner = %d i_loc_m = %d\n",p_m,kk,l2g_block(kk,p,i_loc_m)%p,i_loc_m);
    

    for(int i = 0; i < p_m; i++)
        for(int j = 0; j < n; j++)
        {
            start[i*n+j] = buf[i*n+j];
        }

}

int MPI_Solve(double *a, double *b, double *x,int n,int m,int p,int kk,
    double *buf,double *tmpbuf,double *vecbuf,double *resvec,
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
    resvec=resvec;
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

    int k_bl = n/m, l = n - k_bl*m;
    l=l;
    int is_l = (l==0 ? 0:1);
    int b_loc_rows = get_block_rows(n,m,p,kk);
    int owner;
    int last_owner = (k_bl + is_l - 1)%p;
    DI main_block{},recv_main_block{};
    (void)b_loc_rows;
    main_block=main_block;
    recv_main_block=recv_main_block;
    double eps = 1e-15*N;
    int err{};
    // MPI_Status st;
    eps=eps;
    int rows = get_rows(n,m,p,kk);
    int b_rows = get_block_rows(n,m,p,kk);


    // printf("Printing my loc matrix proc %d:\n",kk);
    // for(int i = 0; i < get_rows(n,m,p,kk); i++)
    // {
    //     for(int j = 0; j < n; j++)
    //     {
    //         printf(" %10.3e",a[i*n+j]);
    //     }
    //     printf("\n");
    // }

    for(int i_glob_m = 0; i_glob_m < k_bl+is_l; i_glob_m++)
    {
        owner = i_glob_m%p;
        int i_loc_m = g2l_block(p,i_glob_m);
        int p_m = (i_glob_m == k_bl ? l:m);
        
        get_block_row(a,buf,n,m,p_m,p,kk,i_loc_m,i_glob_m);
        
        MPI_Bcast(buf,n*p_m,MPI_DOUBLE,owner,MPI_COMM_WORLD);//send block row to all

        // if(kk == main_kk)
        // {
        //     printf("Owner %d printing proc %d buf i_loc_m = %d i_glob_m = %d:\n",owner,kk,i_loc_m,i_glob_m);
        //     for(int i = 0; i < p_m; i++)
        //     {
        //         for(int j = 0; j < n; j++)
        //         {
        //         printf(" %10.3e",buf[i*n+j]);
        //         }
        //         printf("\n");
        //     }
        // }
        
        main_block.norm = 1e64; // сюда будем класть 1/norm
        main_block.num = -1;
        
        if(i_glob_m != k_bl)
        {
            for(int j_loc_m = i_glob_m + kk; j_loc_m < k_bl; j_loc_m += p)
            {
                get_block(buf,block_mm,n,m,0,j_loc_m,i_glob_m);

                // printf("Block[%d,%d] proc %d owner %d in row %d\n",i_glob_m,j_loc_m,kk,owner,i_glob_m);
                // printlxn(block_mm,m,m,m,m);
                if(/*!(inverse1(m,block_mm,invblock_mm,tmpblock_mm))*/inverse(invblock_mm,block_mm,tmpblock_mm,m,eps))
                {
                    N = normofmatrix(invblock_mm,m);
                    if( N < main_block.norm) 
                    {
                        main_block.norm = N;
                        main_block.num = j_loc_m;
                    }
                }
            }// конец проверки строки на вырожденность
        }
        // после реализации прямого хода разкоментить
        else if (kk == last_owner && i_glob_m == k_bl && l!=0)//проверка [l,l] блока
        {
            // printf("Proc %d in ll block\n",kk);
            get_block(buf,block_ll,n,m,0,k_bl,i_glob_m);
            // printf("Block[%d,%d] (l,l) proc %d owner %d in row %d\n",i_glob_m,k_bl,kk,owner,i_glob_m);
            // printlxn(block_ll,l,l,l,l);
            if(/*inverse1(l, block_ll, invblock_ll, tmpblock_ll)*/!inverse(invblock_ll,block_ll,tmpblock_ll,l,eps))
            {
                printf("Block [%d,%d] (block[l,l] in our matrix)  has no inverse after the transformations\n",i_glob_m,k_bl);
                // err = -1;
            }
            else
            {
                main_block.norm = normofmatrix(invblock_ll,l);
                main_block.num = k_bl;
            }
        }//конец проверки [l,l] блока
        MPI_Bcast(&err,1,MPI_INT,last_owner,MPI_COMM_WORLD);
        if(err) return err;

        main_block.norm = 1.0/main_block.norm;

        //обмен процессов, чтобы выяснить у кого главный юлок в строке
        MPI_Allreduce(&main_block,&recv_main_block,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD); // ??
        main_block = recv_main_block;
        // if(kk == main_kk) printf("In proc %d row %d main_block.norm = %10.3e main_block.num = %d\n",kk,i_glob_m, main_block.norm,main_block.num);

        if(kk == main_kk)
        {    
            if((fabs(main_block.norm - 1e64) < eps))
            {   
                
                if(i_glob_m != 0)
                {
                    printf("No inverse matrix in row %d after the transformations in proc %d\n",i_glob_m,owner);
                    // cout<<"\n MATRIX A :\n";
                    // printlxn(a,n,n,n,r);
                    err = -1;
                    
                }
                else
                    printf("No inverse matrix in row %d\n",i_glob_m);

                // clear(block_mm,block_ml,block_ll,tmpblock_mm,tmpblock_ml,tmpblock_ml1,tmpblock_ll,invblock_mm,invblock_ll,diagblock_mm,diaginvblock_mm,vecb_m,vecb_l,tmpvecb_m, tmpvecb_l,colsw);
                err = -1;
            }
        }

        if (main_block.num < 0) return -1;
        // MPI_Bcast(&err,1,MPI_INT,main_kk,MPI_COMM_WORLD);
        // if(err) return err;

        if(main_block.num != i_glob_m && i_glob_m < k_bl)
        {

            // printf("BEFORE SWAP COLUMNS %d %d\n",i_glob_m,main_block.num);
            // print_matrix(a,n,m,p,kk,tmpbuf,n,MPI_COMM_WORLD);
            // printf("\nBUF:\n");
            // printlxn(buf,n,m,rows,n);// мб потом сделать просто вывод buf
            // printf("\n");
            swap_block_columns(a,n,m,rows,i_glob_m,main_block.num);
            swap_block_columns(buf,n,m,p_m,i_glob_m,main_block.num);
            //swap_block_columns(buf,...)
            // printf("AFTER SWAP COLUMNS %d %d\n",i_glob_m,main_block.num);
            // printlxn(a,n,n,n,n);// мб потом сделать просто вывод buf
            // print_matrix(a,n,m,p,kk,tmpbuf,n,MPI_COMM_WORLD);
            // printf("\nBUF:\n");
            // printlxn(buf,n,m,rows,n);
            // printf("\n");
            swap(colsw[i_glob_m],colsw[main_block.num]);
            // if(kk == main_kk) printf("swapped %d %d in row %d\n",i_glob_m,main_block.num,i_glob_m);
        }
        else if(main_block.num == -1 && i_glob_m < k_bl)
        {
            printf("NO mainblock in row %d\n",i_glob_m);
            err = -1;
        }
        // MPI_Bcast(&err,1,MPI_INT,main_kk,MPI_COMM_WORLD);
        if(err) return err;

        //start multiplication

        if(i_glob_m < k_bl )
        {
            get_block(buf,diagblock_mm,n,m,0,i_glob_m,i_glob_m);

            if(/*inverse1(m, diagblock_mm, diaginvblock_mm, tmpblock_mm)*/!(inverse(diaginvblock_mm,diagblock_mm,tmpblock_mm,m,eps)))
            {
                printf("no blocks in row has inverse block i = %d proc %d\n",i_glob_m,kk);
                // printlxn(a,n,n,n,r);
                // CLEAR;
                err = -1;
            }
            // MPI_Bcast(&err,1,MPI_INT,main_kk,MPI_COMM_WORLD);
            // if(err) return err;

            if(kk == owner)
            {
                get_vec_block(b,vecb_m,n,m,i_loc_m,i_glob_m);//надо подумать как переписать чтоб было корректно (ещё не тестил в исх виде)
                mat_x_vector(tmpvecb_m,diaginvblock_mm,vecb_m,m);// double *resvec = mat_x_vector(diaginvblock_mm,vecb_m,m);
                // printf("tmpvecb_m proc %d row %d: \n",kk,i_glob_m);
                // printlxn(tmpvecb_m,m,1,m,m);    
                set_vec_block(b,tmpvecb_m,n,m,i_loc_m,i_glob_m);
                // printf("VECTOR B:\n");
                // for(int i = 0; i < get_rows(n,m,p,kk); i++)
                // {
                //     printf(" %10.3e",b[i]);
                // }
            }

            for(int j_loc_m = i_glob_m + 1; j_loc_m < k_bl; j_loc_m ++)
            {
                // printf("\nJ_LOC_M = %d proc %d i_glob_m = %d\n",j_loc_m,kk,i_glob_m);
                get_block(buf,block_mm,n,m,0,j_loc_m,i_glob_m);

                multiplication(tmpblock_mm,diaginvblock_mm,block_mm,m,m,m);// matmult(tmpblock_mm,diaginvblock_mm,block_mm,m,m,m);// double *resmult = matmult(diaginvblock_mm,block_mm,m,m,m)

                set_block(buf,tmpblock_mm,n,m,0,j_loc_m,i_glob_m);
            }

            if(is_l != 0)
            {
                get_block_ml(buf,block_ml,n,m,l,0);
                multiplication(tmpblock_ml,diaginvblock_mm,block_ml,m,m,l);// matmult(tmpblock_ml,diaginvblock_mm,block_ml,m,m,l);
                set_block_ml(buf,tmpblock_ml,n,m,l,0);
            }

        }
        else if(kk == last_owner && i_glob_m == k_bl)
        {
            get_block(buf,block_ll,n,m,0,i_glob_m,i_glob_m);
            get_vec_block(b,vecb_l,n,m,i_loc_m,i_glob_m);

            if(/*inverse1(l, block_ll, invblock_ll, tmpblock_ll)*/!(inverse(invblock_ll,block_ll,tmpblock_ll,l,eps)))
            {
                printf("ll block has no inverse in proc %d\n",kk);
                // CLEAR; // нужно сделать не выход ретурн -1 а завести флаг по которому потом выйдут одновременно все потоки а то будет DL
                err = -1;
            }
            else
            {
            multiplication(tmpblock_ll,invblock_ll,block_ll,l,l,l);

            mat_x_vector(tmpvecb_l,invblock_ll,vecb_l,l);

            set_block(buf,tmpblock_ll,n,m,0,i_glob_m,i_glob_m);
            set_vec_block(b,tmpvecb_l,n,m,i_loc_m,i_glob_m);
            }

        }
        MPI_Bcast(&err,1,MPI_INT,owner,MPI_COMM_WORLD);
        if(err) return err;


        
        
        
        // for(int j_loc_m = i_glob_m; j_loc_m < k_bl + is_l; j_loc_m++)
        // {
        //     if((j_loc_m - i_glob_m)%p != kk)
        //     {
        //         int col0 = j_loc_m * m;
        //         int size = (j_loc_m == k_bl ? l:m);

        //         for (int ii = 0; ii < p_m; ii++)
        //             memset(buf + ii*n + col0, 0, size * sizeof(double));
        //     }
        // }

        // if(kk != owner)
        // {
        //     printf("\nBEFORE ALLREDUCE Owner %d printing proc %d buf i_loc_m = %d i_glob_m = %d:\n",owner,kk,i_loc_m,i_glob_m);
        //     for(int i = 0; i < p_m; i++)
        //     {
        //         for(int j = 0; j < n; j++)
        //         {
        //         printf(" %10.3e",buf[i*n+j]);
        //         }
        //         printf("\n");
        //     }
        // }

        

        // MPI_Allreduce(buf,tmpbuf,p_m*n,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        // memcpy(buf, tmpbuf, p_m*n*sizeof(double));
        // MPI_Bcast(buf,p_m*n,MPI_DOUBLE,owner,MPI_COMM_WORLD);
        
        set_block_row(a,buf,n,m,p_m,p,kk,i_loc_m,i_glob_m);
        // if(kk == main_kk) printf("\n\nMATRIX A AFTER MULT ROW %d:\n",i_glob_m);
        // print_matrix(a,n,m,p,kk,tmpbuf,n,MPI_COMM_WORLD);
        // if(kk == main_kk) printf("\nVECTOR B AFTER MULT ROW %d:\n",i_glob_m);
        // print_vector(b,n,m,p,kk,vecbuf,n,MPI_COMM_WORLD);
        // if(kk == main_kk) printf("\n");

        if(kk == owner)
            memcpy(resvec,b+i_loc_m*m,p_m*sizeof(double));
        
        MPI_Bcast(resvec,p_m,MPI_DOUBLE,owner,MPI_COMM_WORLD);//правильно ли работае с ll?

        // if(kk == owner)
        // {
        //     printf("MY PART(proc %d i_glob_m = %d) OF VECTOR B: \n",kk,i_glob_m);
        //     for(int i = 0; i < p_m; i++)
        //     {
        //         printf(" %10.3e",resvec[i]);
        //     }
        //     printf("\n");
        // }

        
        

        for(int ii_loc_m = 0; ii_loc_m < b_rows; ii_loc_m++)
        {
            int ii_glob_m = l2g_block(kk,p,ii_loc_m);
            // printf("\nii_glob_m = %d i_glom_m = %d proc %d\n",ii_glob_m,i_glob_m,kk); // правильно

            if(ii_glob_m <= i_glob_m) continue;//если номер локальной строки в глобальной нумерации <= текущего глобального номера то пропускаем
            // printf("\nii_glob_m = %d i_glob_m = %d proc %d\n",ii_glob_m,i_glob_m,kk); // правильно


            if(ii_glob_m != k_bl )
            {
                get_block(a,block_mm,n,m,ii_loc_m,i_glob_m,ii_glob_m); // начал исправлять индексы
                // get_block(a,tmpblock_mm,n,m,ii_loc_m,i_glob_m);
                // printf("\nproc %d take block_mm[%d,%d]\n",kk,ii_glob_m,i_glob_m);
                // printlxn(block_mm,m,m,m,m);

                // memset(tmpblock_mm,0, m*m*sizeof(double));
                // set_block(a,tmpblock_mm,n,m,ii_loc_m,i_glob_m,ii_glob_m);

                get_vec_block(resvec,vecb_m,n,m,0,i_glob_m);//resvec должны переслать часть вектора b из owner всем //вычитание из вектора b block_mm*b
                get_vec_block(b,tmpvecb_m,n,m,ii_loc_m,ii_glob_m);
                vec_mult_sub(tmpvecb_m,block_mm,vecb_m,m);
                set_vec_block(b,tmpvecb_m,n,m,ii_loc_m,ii_glob_m);

                for(int j_loc_m = i_glob_m + 1; j_loc_m < k_bl; j_loc_m++)
                {
                    get_block(buf,invblock_mm,n,m,0,j_loc_m,i_glob_m);
                    get_block(a,diagblock_mm,n,m,ii_loc_m,j_loc_m,ii_glob_m);
                    mat_mult_sub(diagblock_mm,block_mm,invblock_mm,m,m,m);
                    set_block(a,diagblock_mm,n,m,ii_loc_m,j_loc_m,ii_glob_m);
                }

                if (is_l!= 0) 
                {
                    get_block_ml(buf,tmpblock_ml,n,m,l,0);
                    get_block_ml(a,tmpblock_ml1,n,m,l,ii_loc_m);
                    mat_mult_sub(tmpblock_ml1,block_mm,tmpblock_ml,m,l,m);
                    set_block_ml(a,tmpblock_ml1,n,m,l,ii_loc_m);
                }
                // тут вроде все индексы правильные
            }
            else if(kk == last_owner && ii_glob_m == k_bl && l!=0)
            {
                
                get_block_lm(a, block_ml,rows,n, m, l, i_glob_m);


                // get_block_lm(a, tmpblock_ml, rows, m, l, i_loc_m);
                // memset(tmpblock_ml,0,m*l*sizeof(double));
                // set_block_lm(a, tmpblock_ml,rows, n, m, l, i_glob_m);

                // printf("block_lm in col %d(global)  proc %d\n",i_glob_m,kk);
                // printlxn(block_ml,m,l,m,n);
                // printf("\n");

                get_vec_block(resvec,vecb_m,n,m,0,i_glob_m);

                get_vec_block(b,tmpvecb_l,n,m,ii_loc_m,ii_glob_m);
                vec_mult_sub_lm(tmpvecb_l,block_ml,vecb_m,l,m);// vec_mult_sub(tmpvecb_m,block_mm,vecb_m,m);
                set_vec_block(b,tmpvecb_l,n,m,ii_loc_m,ii_glob_m);  // set_vec_block(b,tmpvecb_m,n,m,r);
                
                for(int j_loc_m = i_glob_m + 1; j_loc_m < k_bl; j_loc_m++)
                {
                    get_block(buf,tmpblock_mm,n,m,0,j_loc_m,i_glob_m);
                    get_block_lm(a, tmpblock_ml,rows,n, m, l, j_loc_m);

                    // printf("tmpblock_ml(lxm) in col %d proc %d\n",j_loc_m,kk);
                    // printlxn(tmpblock_ml,m,l,m,m);
                    // printf("\n");
                    
                    mat_mult_sub(tmpblock_ml,block_ml,tmpblock_mm,l,m,m);

                    // printf("AFTER REDUCE tmpblock_ml(lxm) in col %d proc %d\n",j_loc_m,kk);
                    // printlxn(tmpblock_ml,m,l,m,m);
                    // printf("\n");
                    //good


                    // get_vec_block(b,tmpvecb_m,n,m,ii_loc_m,ii_glob_m);//
                    // vec_mult_sub_lm(tmpvecb_m,block_ml,vecb_m,l,m);//
                    set_block_lm(a, tmpblock_ml, rows, n,m, l, j_loc_m);
                    // set_vec_block(b,tmpvecb_m,n,m,r);//set_vec...


                }

                if (is_l != 0)
                {
                    get_block_ml(buf,tmpblock_ml,n,m,l,0);
                    get_block(a,tmpblock_ll,n,m,b_rows-1,k_bl,ii_glob_m);
                    mat_mult_sub(tmpblock_ll,block_ml,tmpblock_ml,l,l,m);
                    set_block(a,tmpblock_ll,n,m,b_rows-1,k_bl,ii_glob_m);
                    //переделать get_block и все его версии
                }
            }

            
        }
        // if(kk == main_kk) printf("\n\nMATRIX A AFTER REDUCE ROW %d:\n",i_glob_m);
        //     print_matrix(a,n,m,p,kk,tmpbuf,n,MPI_COMM_WORLD);
        //     if(kk == main_kk) printf("\nVECTOR B AFTER REDUCE ROW %d:\n",i_glob_m);
        //     print_vector(b,n,m,p,kk,vecbuf,n,MPI_COMM_WORLD);
        //     if(kk == main_kk) printf("\n");//что то не так в прямом ходе посмотреть чтательно индексы
        

    }//end straight algo(вроде всё работает)

    //start reverse algo
    memset(tmpbuf,0,n*sizeof(double));

    for(int i_glob_m = k_bl + is_l - 1; i_glob_m >= 0; i_glob_m--)
    {
        int i_loc_m = g2l_block(p,i_glob_m);
        owner = i_glob_m%p;
        int p_m = (i_glob_m == k_bl ? l:m);

        if(kk == owner)
        {
            for(int s = 0; s < p_m; s++)
            {
                if(p_m != m) vecb_l[s] = b[i_loc_m*m + s];
                else vecb_m[s] = b[i_loc_m*m + s];
            }

            for(int j_loc_m = i_glob_m + 1; j_loc_m < k_bl + is_l ; j_loc_m++)
            {
                int col_size = (j_loc_m == k_bl ? l:m);

                for(int ii = 0; ii < p_m; ii++)
                {
                    double sum{};
                    for(int jj = 0; jj < col_size; jj++)
                    {
                        sum += a[ (i_loc_m*m + ii)*n + j_loc_m*m + jj] * tmpbuf[j_loc_m*m + jj];
                    }

                    if(p_m != m) vecb_l[ii] -= sum;
                    else vecb_m[ii] -= sum;

                }

            }
        }

        if(p_m != m) MPI_Bcast(vecb_l,p_m,MPI_DOUBLE,owner,MPI_COMM_WORLD);
        else MPI_Bcast(vecb_m,p_m,MPI_DOUBLE,owner,MPI_COMM_WORLD);

        if(p_m != m) for(int t = 0; t < l; t++) tmpbuf[t+i_glob_m*m] = vecb_l[t];
        else for(int t = 0; t < m; t++) tmpbuf[t+i_glob_m*m] = vecb_m[t];
    }

    //получили вектор х в первой строке tmpbuf

    //начинаем перестановку

    // printf("colsw in proc %d:\n",kk);
    for(int i = 0; i < k_bl; i++)
    {
        // printf(" %d",colsw[i]);
        if (colsw[i] != i )
        {
            swap_block_columns(tmpbuf,n,m,1,i,colsw[i]);
            swap(colsw[i],colsw[colsw[i]]);
        }
    }
    // printf("\n");

    for(int ii_loc_m = 0; ii_loc_m < rows; ii_loc_m++)
    {
        int ii_glob_m = l2g(m,kk,p,ii_loc_m);

        x[ii_loc_m] = tmpbuf[ii_glob_m];
    }
    
    

    return 0;
}