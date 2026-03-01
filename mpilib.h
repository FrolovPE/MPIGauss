#ifndef MPILIB_H
#define MPILIB_H

#include <iostream>
#include <chrono>
#include <unistd.h>
#include <cstring>
#include <vector>
#include <pthread.h>
#include <sched.h>
#include <sys/sysinfo.h>
#include "mytime.h"
#include "mpi.h"




#define EPS64 1e-64
#define CLEAR     clear(block_mm,block_ml,block_ll,tmpblock_mm,tmpblock_ml,tmpblock_ml1,tmpblock_ll,invblock_mm,invblock_ll,diagblock_mm,diaginvblock_mm,vecb_m,vecb_l,tmpvecb_m, tmpvecb_l,colsw)
using namespace std;

struct DI
{
    double norm;
    int num;
};

void printlxn(const double *a, int size, int l, int n, int r);
int readarray(double *a, int n, const char* filename);
int read_array(FILE *fp,double *buf,int &printed,int n, double totalsize);
double f (int s , int n , int i , int j);
void init(double *a, double (*f)(int,int,int,int), int n, int s);
int solution(int n, int m, double *a, double *b, double *x,
    double *block_mm, double *block_ml, double *block_ll,double *invblock_mm, double *diaginvblock_mm, 
    double *invblock_ll,double *diagblock_mm,
    int *colsw,double *vecb_m,double *vecb_l,
    double *tmpblock_mm,double *tmpblock_ml,double *tmpblock_ml1,double *tmpblock_ll,double *tmpvecb_m,double *tmpvecb_l,const double eps);
double vectornorm(double *a , int n);
void mpi_norm_matrix(double *a,double &N,int n,int m, int p, int k,MPI_Comm com);
void mat_x_vector(double *res,double *a, double *b, int n);
void vectorsub(double *res,double *a,double *b, int n);
void residuals(double &r1,double &r2,double *a,double *b,double *x,double *realx,int n,double *Ax, double *Ax_b, double *x_realx);
void mpi_residuals(double &r1,double &r2,double *a,double *b,double *x,double *realx,int n,int m,int p, int k, double *Ax, double *Ax_b, double *x_realx,MPI_Comm com);
void report(char *title, int task, double r1, double r2 ,double t1,  double t2 ,int s, int n , int m,int p );
void matmult(double *res,double *a, double *b, int n, int m,int l);
void get_block(double *a, double *b, int n, int m, int i, int j);
void set_block(double *a, double *b, int n, int m, int i, int j);
double normofmatrix(double *a , int size);
double* inverse(double *result,const double* A, int size,double eps);
double* inverse1(const double* a_in, int n,const double eps );
void blocksize(int i, int j,int n,int m, int &r, int &h);
void swap_block_columns(double *a, int n,int m, int i, int j);
void swap_block_vec(double *a, int n,int m, int i, int j);
void get_vec_block(double *b,double *block,int n, int m, int i);
void mat_mult_sub(double* Result, double* Block_A, double* Block_B, const int row_A, const int col_B, const int col_row);
void set_vec_block(double *b,double *block,int n, int m, int i);
void get_block_ml(double *a, double *b, int n, int m,int l, int i);
void set_block_ml(double *a, double *b, int n, int m,int l, int i);
void get_block_lm(double *a, double *b, int n, int m,int l, int j);
void set_block_lm(double *a, double *b, int n, int m,int l, int j);
void matsub(double *res,double *a, double *b, int m,int l);
void vec_eq(double *x, double *b, int m);
void vec_sub(double *x, double *b, int m);
void vec_mult_sub(double* Result, double* A, double* vec, int m);
void vec_mult_sub_lm(double* Result, double* A, double* vec, int l,int m);
void multiplication(double* Result, double* Block_A, double* Block_B, const int row_A, const int col_row,const int col_B);
// void* parallelSolve(void* ptr);
// void* parallelSolve1(void* ptr);
// void pllinit_matrix(double *a, int s,int n , int m , int k, int p);
// void pllinit_vectorb(double *b,double *a,int n , int m , int k, int p);
// void clear(double *block_mm,double *block_ml,double *block_ll,double *tmpblock_mm,double *tmpblock_ml,double *tmpblock_ml1,double *tmpblock_ll,double *invblock_mm,double *invblock_ll,double *diagblock_mm,double *diaginvblock_mm,double *vecb_m,double *vecb_l,double *tmpvecb_m, double *tmpvecb_l,int *colsw);
int l2g(/*int n,*/ int m, int k, int p, int i_loc);
int g2l(/*int n,*/ int m, /*int k,*/ int p, int i_glob);
int l2g_block(/*int n, int m,*/ int k, int p, int i_loc_m);
int g2l_block(/*int n, int m, int k,*/ int p, int i_glob_m);
int get_k(int n, int m, int p, int i_glob);
int get_block_rows(int n, int m, int p, int k);
int get_max_block_rows(int n, int m, int p);
int get_rows(int n, int m, int p, int k);
void init_matrix(double *a, int n, int m, int p, int k, double (*f)(int,int,int,int),int s);
void init_vector_b(double *a,double *b,int n, int m, int p, int k);
int read_matrix(double *a, int n, int m, int p, int k, const char *name,double *buf, MPI_Comm com);
void print_matrix(double *a, int n, int m, int p, int k, double *buf, int max_print, MPI_Comm com);
void print_vector(double *b,int n,int m, int p, int k,double *buf,int max_print,MPI_Comm com);
int print_array(double *a, int n, int m, int printed_rows,int max_print);
void matrix_mult_vector(double *a, double *bvec,double *c, int n, int m, int k, int p,MPI_Comm com);
int MPI_Solve(double *a, double *b, double *x,int n,int m,int p,int kk,double *buf,double *vecbuf,double *block_mm, double *block_ml, double *block_ll,double *invblock_mm, double *diaginvblock_mm, 
    double *invblock_ll,double *diagblock_mm,
    int *colsw,double *vecb_m,double *vecb_l,
    double *tmpblock_mm,double *tmpblock_ml,double *tmpblock_ml1,double *tmpblock_ll,double *tmpvecb_m,double *tmpvecb_l,MPI_Comm com);
void get_block_row(double *a,double *buf,int n, int m, int p,int kk,int i_loc_m);
#endif