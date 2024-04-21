// Search for a good implementation of a ratational an LBox
// Uses greedy decomposition

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

// Number of 32-bit words
#define WIDTH 4

// Clang extension for fixed-size integers
typedef uint32_t word32;
typedef unsigned _BitInt(WIDTH*32) word;

#define MOD(n) ((((n)%(WIDTH*32))+WIDTH*32)%(WIDTH*32))
#define rot(w,n) ((w<<MOD(n))|(w>>(32*WIDTH-MOD(n))))

word rand_word() {
  word x = 0;
  for (int i=0; i<WIDTH; i++) {
    x<<=16;
    x += rand();
    x<<=16;
    x += rand();
  }

  return x;
}

int weight (word x) {
  int w = 0;
  for (int i=WIDTH-1; i>=0; i--) {
    w += __builtin_popcount((uint32_t)(x>>(32*i)));
  }

  return w;
}

int naive(word x[], int n) {
  int w = 0;
  for (int i=0; i<n; i++) {
    w += x[i]?weight(x[i])-1:0;
  }
  return w;
}

int first_bit (word x) {
  word m = 1;
  for (int i=0; i<32*WIDTH; i++) {
    if (x & m)
      return i;
    m <<= 1;
  }

  return -1;
}

void print_word(word x) {
  for (int i=WIDTH-1; i>=0; i--) {
    printf ("%08x", (uint32_t)(x>>(32*i)));
  }
  printf("(%2i)", weight(x));
}

void print_vector(word x[], int n) {
  printf("[");
  for (int i=0; i<n; i++) {
    printf(" ");
    print_word(x[i]);
  }
  printf(" ]");
}

struct decomp {
  word mask;
  word rem;
  int rot;
};

// Decompose word x as rem ^ mask ^ mask<<<rot 
struct decomp decompose (word x) {
  struct decomp ret;
  int best_w = 32*WIDTH;
  
  for (int i=1; i<32*WIDTH; i++) {
    word mask = 0;
    word tmp = x ^ mask ^ rot(mask,i);
    word z;
    while ((z = (tmp & rot(tmp,-i)))) {
      if (z & ~rot(z,i)) {
	mask |= z & ~rot(z,i);
      } else {
	mask |= z & ~(z-1);
      }

      /* assert(weight(mask & ~x) == 0); */
      /* assert(weight(rot(mask,i) & ~x) == 0); */
      /* assert(weight(rot(mask,i) & mask) == 0); */
      
      tmp = x ^ mask ^ rot(mask,i);
    }

    int w = weight(mask)+weight(x ^ mask ^ rot(mask,i));
    
    if (w < best_w) {
      best_w = w;
      ret.mask = mask;
      ret.rem  = x ^ mask ^ rot(mask,i);
      ret.rot  = i;
    }
  }

  return ret;
}

// Try several ways to decompose a set of values
// Uses recursivity to decompose a larger set
// Naive evaluation of each node using weight
int recursive_decompose_greedy(word x0[], int n) {
  word x[n];
  memcpy(x, x0, sizeof(x));
  word y[n+1];
  int y_n = n;
  int cost = 0;
  int w0 = naive(x, n);
  int best = w0;

  // x[i] = 1 <<< s
  for (int i=0; i<n; i++) {
    if (weight(x[i]) == 1) {
      memcpy(y,       x,       i*sizeof(x[0]));
      memcpy(y+i, x+i+1, (n-i-1)*sizeof(x[0]));
      memcpy(x, y, sizeof(y));
      n = n-1;
    }
  }

  // x[i] = x[j] <<< s
  for (int i=0; i<n; i++) {
    for (int j=i+1; j<n; j++) {
      for (int s=0; s<32*WIDTH; s++) {
	if (x[i] == rot(x[j],s)) {
	  memcpy(y,       x,       i*sizeof(x[0]));
	  memcpy(y+i, x+i+1, (n-i-1)*sizeof(x[0]));
	  memcpy(x, y, sizeof(y));
	  n = n-1;
	  w0 = naive(x, n);

	  j = n;
	  break;
	}
      }
    }
  }
  
  // x[i] = rem ^ mask ^ rot(mask,s)
  for (int i=0; i<n; i++) {
    if (weight(x[i]) > 3) {
      struct decomp d = decompose(x[i]);
      int w = w0-(weight(x[i])-1)+(d.rem?weight(d.rem)-1:0)+(weight(d.mask)-1)+1+(weight(d.rem)>0);
      if (w < best) {
	memcpy(y, x, n*sizeof(x[0]));
	y[i] = d.mask;
	y[n] = d.rem;
	y_n = d.rem? n+1: n;
	best = w;
	cost = 1+(d.rem!=0);
      }
    }
  }

  // x[i] ^= x[j] <<< s
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      if (i==j)
	continue;
      for (int s=0; s<32*WIDTH; s++) {
  	int w = w0+1-weight(x[i])+weight(x[i]^rot(x[j],s));
  	if (w < best) {
  	  memcpy(y, x, n*sizeof(x[0]));
  	  y[i] ^= rot(x[j],s);
  	  y_n = n;
  	  best = w;
	  cost = 1;
  	}
      }
    }
  }

  // x[i] ^= mask; x[j] ^= rot(mask,s)
  for (int i=0; i<n; i++) {
    for (int j=i+1; j<n; j++) {
      for (int s=0; s<32*WIDTH; s++) {
  	word m = x[i] & rot(x[j],s);
  	int w = w0+1-weight(m); // false if x[i] == mask or x[j] == mask
  	if (w < best) {
  	  memcpy(y, x, n*sizeof(x[0]));
  	  y[i] ^= m;
  	  y[j] ^= rot(m,-s);
  	  y[n] = m;
  	  y_n = n+1;
  	  best = w;
	  cost = 2;
  	}
      }
    }
  }

  int w = w0;
  if (best < w0)
    w = recursive_decompose_greedy(y, y_n);
  else
    cost = 0;
  
  return w+cost;
}


// Try several ways to decompose a set of values
// Uses recursivity to decompose a larger set
// Advanced evaluation of each node using recursive_decompose_greedy
int recursive_decompose(word x[], int n) {
  word y[n+1];
  int y_n = n;
  word tmp_y[n+1];
  int tmp_yn = n;
  int cost = 0;
  int w0 = naive(x, n);
  int best = w0;

  char BUF[128] = {0};
  
  // x[i] = 1 <<< s
  for (int i=0; i<n; i++) {
    if (weight(x[i]) == 1) {
      memcpy(y,       x,       i*sizeof(x[0]));
      memcpy(y+i, x+i+1, (n-i-1)*sizeof(x[0]));
      int ret = recursive_decompose(y, n-1);
      //print_vector(x,n);
      //printf(" [R%i = 1<<<s]\n", i);
      return ret;
    }
  }

  // x[i] = x[j] <<< s
  for (int i=0; i<n; i++) {
    for (int j=i+1; j<n; j++) {
      for (int s=0; s<32*WIDTH; s++) {
	if (x[i] == rot(x[j],s)) {
	  memcpy(y,       x,       i*sizeof(x[0]));
	  memcpy(y+i, x+i+1, (n-i-1)*sizeof(x[0]));
	  int ret = recursive_decompose(y, n-1);
	  //print_vector(x,n);
	  //printf(" [R%i = R%i<<<%i]\n", i, j, s);
	  return ret;
	}
      }
    }
  }
  
  // x[i] = rem ^ mask ^ rot(mask,s)
  for (int i=0; i<n; i++) {
    if (weight(x[i]) > 3) {
      // struct decomp d = decompose(x[i]);
      // Unroll function to keep several decompositions
    
      for (int r=1; r<32*WIDTH; r++) {
	word mask = 0;
	word tmp = x[i] ^ mask ^ rot(mask,r);
	word z;
	while ((z = (tmp & rot(tmp,-r)))) {
	  if (z & ~rot(z,r)) {
	    mask |= z & ~rot(z,r);
	  } else {
	    mask |= z & ~(z-1);
	  }

	  tmp = x[i] ^ mask ^ rot(mask,r);
	}

	if (weight(mask) > 1) {
	  word rem = x[i] ^ mask ^ rot(mask,r);
	  memcpy(tmp_y, x, n*sizeof(x[0]));
	  tmp_y[i] = mask;
	  tmp_y[n] = rem;
	  tmp_yn = rem? n+1: n;

	  int wr = recursive_decompose_greedy(tmp_y, tmp_yn);
	  if (wr+1+(rem!=0)<best) {
	    cost = 1+(rem!=0);
	    best = wr+cost;
	    memcpy(y, tmp_y, sizeof(y));
	    y_n = tmp_yn;
	    if (rem) {
	      snprintf(BUF, sizeof(BUF), "R%i = R%i ^ R%i ^ R%i<<<%i",
		       i, n, i, i, r);
	    } else {
	      snprintf(BUF, sizeof(BUF), "R%i = R%i ^ R%i<<<%i",
		       i, i, i, r);
	    }
	  }
	}
      }
    }
  }

  // x[i] ^= x[j] <<< s
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      if (i==j)
	continue;
      for (int s=0; s<32*WIDTH; s++) {
	memcpy(tmp_y, x, n*sizeof(x[0]));
	tmp_y[i] ^= rot(x[j],s);
	tmp_yn = n;
	int wr = recursive_decompose_greedy(tmp_y, tmp_yn);
	  
  	if (wr+1 < best) {
	  cost = 1;
	  best = wr+cost;
	  memcpy(y, tmp_y, sizeof(y));
	  y_n = tmp_yn;
	  snprintf(BUF, sizeof(BUF), "R%i ^= R%i<<<%i",
		   i, j, s);
  	}
      }
    }
  }

  // x[i] ^= mask; x[j] ^= rot(mask,s)
  for (int i=0; i<n; i++) {
    for (int j=i+1; j<n; j++) {
      for (int s=0; s<32*WIDTH; s++) {
  	word m = x[i] & rot(x[j],s);
	if (weight(m) > 1) {
	  memcpy(tmp_y, x, n*sizeof(x[0]));
	  tmp_y[i] ^= m;
	  tmp_y[j] ^= rot(m,-s);
	  tmp_y[n] = m;
	  tmp_yn = n+1;
	  int wr = recursive_decompose_greedy(tmp_y, tmp_yn);

	  if (wr+2 < best) {
	    cost = 2;
	    best = wr+cost;
	    memcpy(y, tmp_y, sizeof(y));
	    y_n = tmp_yn;
	    snprintf(BUF, sizeof(BUF), "R%i ^= R%i; R%i ^= R%i<<<%i",
		     i, n, j, n, s);
	  }
	}
      }
    }
  }

   // x[i] ^= 1 <<< s (random noise) 
   for (int i=0; i<n; i++) {
     for (int s=0; s<32*WIDTH; s++) { 
       memcpy(tmp_y, x, n*sizeof(x[0])); 
       tmp_y[i] ^= ((word)1)<<s; 
       tmp_yn = n; 
       int wr = recursive_decompose_greedy(tmp_y, tmp_yn); 
	  
       if (wr+1 < best) { 
   	cost = 1; 
   	best = wr+cost; 
   	memcpy(y, tmp_y, sizeof(y)); 
   	y_n = tmp_yn; 
   	snprintf(BUF, sizeof(BUF), "R%i ^= 1<<<%i", 
   		 i, s); 
       } 
     } 
   } 

   // x[i] ^= 1 <<< s ^ 1 <<< t (random noise) 
   for (int i=0; i<n; i++) { 
     for (int s=0; s<32*WIDTH; s++) { 
       for (int t=s+1; t<32*WIDTH; t++) { 
   	memcpy(tmp_y, x, n*sizeof(x[0])); 
   	tmp_y[i] ^= ((word)1)<<s; 
   	tmp_y[i] ^= ((word)1)<<t; 
   	tmp_yn = n; 
   	int wr = recursive_decompose_greedy(tmp_y, tmp_yn); 
	  
   	if (wr+2 < best) { 
   	  cost = 2; 
   	  best = wr+cost; 
   	  memcpy(y, tmp_y, sizeof(y)); 
   	  y_n = tmp_yn; 
   	  snprintf(BUF, sizeof(BUF), "R%i ^= 1<<<%i ^ 1<<<%i", 
   		   i, s, t); 
   	} 
       } 
     } 
   } 
  
  printf("[%2i/%2i/%2i]  ", best, recursive_decompose_greedy(x, n), w0);
  print_vector(x, n);
  printf("\n");
  /* assert(best <= recursive_decompose_greedy(x, n)); */
  int w = w0;
  if (best < w0)
    w = recursive_decompose(y, y_n);
  else
    cost = 0;
  
  print_vector(x,n);
  printf(" (%i/%i) [%s]\n", w+cost, best, BUF);
  return w+cost;
}



int main (int argc, char *argv[]) {
    

    if (argc > 1)
    {



      word x[1];
    
      if (strlen(argv[1]) != WIDTH*8) {
      BAD:
	printf ("Unable to parse input: '%s'\n", argv[1]);
	return -1;
      }
      x[0] = 0;
      for (int j=0; j<WIDTH*8; j++) {
	char c = argv[1][j];
	if (c >= 'a')
	  c = c-'a'+'A';
	if (c >= '0' && c <= '9')
	  x[0] = (x[0]<<4) + (c-'0');
	else if (c >= 'A' && c <= 'F')
	  x[0] = (x[0]<<4) + (c-'A'+10);
	else
	  goto BAD;

    }


    
    int w = recursive_decompose(x, 1);
    printf("Final cost: %i\n\n", w);
}
    else {
    // Test with random matrices
    srand(42);
    word r = rand_word();

    
      printf("Seaching implementation of ");
      print_word(r);
      printf("\n");
    
      int w = recursive_decompose(&r, 1);
      printf("Final cost: %i\n\n", w);
     
    
  }
    
}
