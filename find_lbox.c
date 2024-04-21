#define _XOPEN_SOURCE 500

#include <stdio.h>
#include <stdlib.h>
#include <x86intrin.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <inttypes.h>

#include <sys/types.h>
#include <unistd.h>
#include <time.h>
#include <wmmintrin.h>

typedef uint64_t word;
typedef unsigned __int128 dword;
typedef uint32_t half_word;
#define POPCOUNT __builtin_popcountll

// Parameters for the search
// See examples in the README file

#define WORD_SIZE 32
#define WORD_NB 4
#define TARGET_WEIGHT 20
#define ITER 10
// Number of XOR per 32-bit words
#define MAX_XOR 6
// #define COMBINE3ROWS



// Internal macros

#define MATRIX_SIZE (WORD_SIZE*WORD_NB)
#define NB_MASK ((1U<<WORD_NB)-1)
#if MATRIX_SIZE < 128
  #define TRUNC(x) ((x)&((((dword)1)<<MATRIX_SIZE)-1))
#else
  #define TRUNC(x) (x)
#endif

/* Binary matrix */

typedef struct {
  unsigned int m;
  unsigned int n;
  unsigned int bn;
  dword** data;
} matrix;

int gauss_reduce (matrix *M, unsigned int cols, int last);

dword rot (dword x, int r) {
  r = r % MATRIX_SIZE;
  if (r == 0)
    {
    return x;
    }
  x = ((x << r) | (x >> (MATRIX_SIZE - r)));
  return TRUNC(x);
}

void print_row(dword x) {
  for (int i=0; i<MATRIX_SIZE/8; i++) {
    printf ("%02x", (int)((x>>(MATRIX_SIZE-8*i-8))&0xff));
  }
}

int w (dword x) {
  dword xx = 0;
  for (int i=0; i<WORD_NB; i++)
    xx |= x>>i;
  dword mm = 0;
  for (int i=0; i<WORD_SIZE; i++)
    mm = (mm<<WORD_NB) | 1;
  
  xx &= mm;
  
  word xx0 = (word) xx;
  word xx1 = (word) (xx >> 64);
  return POPCOUNT (xx0) + POPCOUNT (xx1);
}

dword invert (dword x, matrix M) {
  if (!x)
    return 0;

  for (int i = 0; i < MATRIX_SIZE; i++) {
    M.data[i][0] = rot (x, i);
    M.data[i][1] = ((dword)1) << i;
  }

  int g = gauss_reduce (&M, MATRIX_SIZE, 0);
  if (!g) {
    return 0;
  }

  return M.data[0][1];
}

void Permute_IS(dword M, int p[], matrix *M_IS) {
#define GET_M(i,j) (( ((i)<WORD_SIZE? rot(M,j): (dword)1<<j) >>(WORD_NB*(i%WORD_SIZE))) & NB_MASK )
  for (int j=0; j<WORD_SIZE*WORD_NB; j++) {
    M_IS->data[j][0] = M_IS->data[j][1] = 0;
    for (int i=0; i<2*WORD_SIZE; i++) {
      M_IS->data[j][i/WORD_SIZE] |= GET_M(p[i],j) << (WORD_NB*(i%WORD_SIZE));
    }
  }
#undef GET_M
}

void random_perm(int p[]) {
  for (int i=0; i<2*WORD_SIZE; i++) {
    int j = rand() % (i+1);
    p[i] = p[j];
    p[j] = i;
  }
}

// Parameters:
// * M:   input matrix
// * Ma:  scratch space for a matrix
// * min: minumum branch number
// * it:  number of iterations (log2)
int Prange (dword M, matrix *Ma, int min, int it)
{
  int p_min = MATRIX_SIZE;
  
  for (int t = 0; t < (1 << it); t++)
    { 
      int p[2*WORD_SIZE];
      random_perm(p);
      Permute_IS(M, p, Ma);

    gauss_reduce(Ma, MATRIX_SIZE, 0);


     /* 1 row */
    for (int i = 0; i < WORD_SIZE; i++) 
      {
      for (int j = 1; j < 1<<WORD_NB; j++)
        {
        dword alpha_g = 0;
        dword alpha_d = 0;
        for (int dd = 0; dd < WORD_NB; dd++)
          {
          if ((j >> dd) % 2 == 1)
            {
            alpha_g ^= Ma->data[WORD_NB * i + dd][1];
            alpha_d ^= Ma->data[WORD_NB * i + dd][0];
            }
          }
        if (w (alpha_g) + w (alpha_d) < p_min)
      	  {
      	  p_min =  w (alpha_g) + w (alpha_d);
      	  }
        }
      }
      
    if (p_min <= min)
      {
      return p_min;
      }
      
    /* 2 rows */
    for (int i = 0; i < WORD_SIZE-1; i++) 
      {
      for (int j = 1; j < 1<<WORD_NB; j++)
        {
        dword alpha_g = 0;
        dword alpha_d = 0;
        for (int dd = 0; dd < WORD_NB; dd++)
          {
          if ((j >> dd) % 2 == 1)
            {
            alpha_g ^= Ma->data[WORD_NB * i + dd][1];
            alpha_d ^= Ma->data[WORD_NB * i + dd][0];
            }
          }
          
          for (int i2=i + 1; i2<WORD_SIZE; i2++) 
            {
            for (int j2 = 1; j2 < 1<<WORD_NB; j2++)
              {
              for (int dd2 = 0; dd2 < WORD_NB; dd2++)
                {
                if ((j2 >> dd2) % 2 == 1)
                  {
                  alpha_g ^= Ma->data[WORD_NB * i2 + dd2][1];
                  alpha_d ^= Ma->data[WORD_NB * i2 + dd2][0];
                  }
                }
                    
              if (w (alpha_g) + w (alpha_d) < p_min)
      	        {
      	        p_min =  w (alpha_g) + w (alpha_d);
      	        }
      	      }
      	    }
      	  }  
        }
      
        /* 3 rows */
#ifdef COMBINE3ROWS
    for (int i=0; i<WORD_SIZE-2; i++) 
      {
      for (int j = 1; j < 1<<WORD_NB; j++)
        {
        dword alpha_g = 0;
        dword alpha_d = 0;
        for (int dd = 0; dd < WORD_NB; dd++)
          {
          if ((j >> dd) % 2 == 1)
            {
            alpha_g ^= Ma->data[WORD_NB * i + dd][1];
            alpha_d ^= Ma->data[WORD_NB * i + dd][0];
            }
          }
          
          for (int i2=i + 1; i2<WORD_SIZE-1; i2++) 
            {
            for (int j2 = 1; j2 < 1<<WORD_NB; j2++)
              {
              for (int dd2 = 0; dd2 < WORD_NB; dd2++)
                {
                if ((j2 >> dd2) % 2 == 1)
                  {
                  alpha_g ^= Ma->data[WORD_NB * i2 + dd2][1];
                  alpha_d ^= Ma->data[WORD_NB * i2 + dd2][0];
                  }
                }
                for (int i3=i2 + 1; i3<WORD_SIZE; i3++) 
                  {
                  for (int j3 = 1; j3 < 1<<WORD_NB; j3++)
                    {
                    for (int dd3 = 0; dd3 < WORD_NB; dd3++)
                      {
                      if ((j3 >> dd3) % 2 == 1)
                        {
                        alpha_g ^= Ma->data[WORD_NB * i3 + dd3][1];
                        alpha_d ^= Ma->data[WORD_NB * i3 + dd3][0];
                        }
                      } 
                    
                    if (w (alpha_g) + w (alpha_d) < p_min)
      	              {
      	              p_min =  w (alpha_g) + w (alpha_d);
      	              }
      	            }
      	          }
      	        }
      	      }
      	    }
      	  }
#endif
    
  if (p_min <= min)
    {
    return p_min;
    }
  }
  //printf ("p = %d\n", p_min);
  return p_min;
}

int main(int argc, char* argv[]) 
{
  int do_random = 0;
  dword do_one = 0;

  if (argc > 2) {
    goto USAGE;
  } else if (argc == 2) {
    if (strcmp(argv[1], "-r") == 0) {
      do_random = 1;
    } else {
      if (strlen(argv[1]) != MATRIX_SIZE/4) {
      USAGE:
	printf ("Usage: %s [-r] [mat]\n", argv[0]);
	return -1;	
      }
      for (int j=0; j<MATRIX_SIZE/4; j++) {
	char c = argv[1][j];
	if (c >= 'a')
	  c = c-'a'+'A';
	if (c >= '0' && c <= '9')
	  do_one = (do_one<<4) + (c-'0');
	else if (c >= 'A' && c <= 'F')
	  do_one = (do_one<<4) + (c-'A'+10);
	else
	  goto USAGE;
      }
    }
  }

  int seed = 1337 * getpid() + time(NULL);
  srand(seed);

#pragma omp parallel
  {
  dword  Mdata[2*MATRIX_SIZE];
  dword* Mptr[MATRIX_SIZE];
  matrix Ma = {.m = MATRIX_SIZE, .n = 2*MATRIX_SIZE, .bn = 2,
		.data = Mptr };
  for (int i = 0; i < MATRIX_SIZE; i++)
    Ma.data[i] = Mdata + 2 * i;
      
          
          
  for (unsigned long long test = 0; test < (do_one? 1: (1 << 20)); test++)
    {
      dword M;
      dword x[MAX_XOR] = {1};
      int r[MAX_XOR] = {0};

      if (do_one) {
	// Matrix given as argument
	M = do_one;
      }
      else if (do_random)
	{
	  // Random test
	  M = (((dword) random ()) << 96) ^ (((dword) random ()) << 64) ^ (((dword) random ()) << 32) ^ ((dword) random ());
	  M = TRUNC(M);
	}
      else
	{
	// Normal test (candidates generated by combinations of XOR and rotations)
	for (int n = 1; n < MAX_XOR; n++) 
	  {
	    // Bias towards recent values
	    int i = n - 1 - __builtin_ctz (random ());
	    if (i < 0) i = 0;
	    int j = n - 1 - __builtin_ctz ( random());
	    if (j < 0) j = 0;
	    r[n] = random () % MATRIX_SIZE;
	    x[n] = x[i] ^ rot (x[j], r[n]);
	  }
	M = x[MAX_XOR-1];
      }

    dword I = invert(M, Ma);
    if (!I)
#pragma omp critical
      {
      test--;
      if (do_one) {
	printf ("Not invertible\n");
      }
      }
    
    else 
      {
      int p_min = Prange (M, &Ma, TARGET_WEIGHT, ITER);
      if (do_one || p_min > TARGET_WEIGHT)
        {
        printf("M = ");
	print_row(M);
        printf(" I = ");
	print_row(I);
        printf ("\nminimum weight found = %d\n", p_min);
	if (!do_random && !do_one) {
	  printf("Operations: ");
	  for (int i=0; i<MAX_XOR; i++) {
	    printf (" (%2i) [", r[i]);
	    print_row(x[i]);
	    printf ("]");
	  }
	  printf ("\n\n\n");
	}
        fflush(stdout);
        }
      }
    }
  }
  return 0;
}


//// Gaussian pivot to invert matrix ////

int gauss_reduce (matrix *M, unsigned int cols, int last) {
  /* column index */
  dword mask = 0;
  unsigned int bk = 0;
  unsigned int wk = -1;
  /* M[i][wk] & mask is M_{i,k} */
  /* mask is 1 << bk */

  /* row index */
  unsigned int l = 0;

  for (unsigned int k = 0; k < cols && l < M->m; k++) 
    {
    mask <<= 1;
    
    bk++;
    if ((!mask)) 
      {
      wk++;
      bk = 0;
      mask = 1;
      }
    unsigned int pos = (unsigned int) -1;
    for (unsigned int i = l; i < M->m - last; i++) 
      {
      if (M->data[i][wk] & mask) 
        {
        pos = i;
        break;
        }
      }

    if (pos != (unsigned int) - 1) 
      {
      if (l != pos) 
        {
	dword *t = M->data[pos];
	M->data[pos] = M->data[l];
	M->data[l] = t;
        }

      dword *y = M->data[l];

      for (unsigned int i = 0; i < M->m; i++) 
        {
        if (i == l)
          continue;

        if (M->data[i][wk] & mask) 
          {
          dword *x = M->data[i];

          for (unsigned int j = wk; j < M->bn; j++)
            x[j] ^= y[j];
          } 
        }

      l++; 
    }
  }

  return (l == M->m);
}


