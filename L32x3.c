#include <stdio.h>
#include <stdlib.h>

#include <stdint.h>

#include <inttypes.h>

typedef uint32_t word;

word rot (word x, int r) {
  r = r % 32;
  if (r == 0)
    {
    return x;
    }
  return ((x << r) | (x >> (32 - r)));
}


// Implementation of L_32x3 LBox
void M (word *x, word *y, word *z)
{
  word a, b, c, d, e, f;
  
  a = *x ^ rot (*y, 17);
  b = *y ^ rot (*z, 17);
  c = *z ^ rot (*x, 18);
  
  d = a ^ rot (c, 30);
  e = b ^ rot (a, 31);
  c = c ^ rot (b, 31);
  
  a = d ^ rot (*x, 21);
  b = e ^ rot (*y, 21);
  f = c ^ rot (*z, 21);
  
  d = d ^ rot (b, 30);
  e = e ^ rot (f, 30);
  c = c ^ rot (a, 31);
  
  a = d ^ rot (c, 23);
  b = e ^ rot (d, 24);
  f = c ^ rot (e, 24);
  
  a = a ^ rot (d, 26);
  b = b ^ rot (e, 26);
  c = f ^ rot (c, 26);
  
  *x = a;
  *y = b;
  *z = c;
}

// Implementation of L_32x3 LBox inverse
void I (word *x, word *y, word *z)
{
  word a, b, c, d, e, f, g, h, i, j, k, l;

  a = rot (*y, 14) ^ rot (*y, 26);
  b = rot (*z, 14) ^ rot (*z, 26);
  c = rot (*x, 15) ^ rot (*x, 27);
  
  d = rot (a, 8);
  e = rot (b, 8);
  f = rot (c, 8);
  
  a = a ^ rot (a, 1);
  b = b ^ rot (b, 1);
  c = c ^ rot (c, 1);
  
  g = rot (*z, 6);
  h = rot (*x, 7);
  i = rot (*y, 7);
  
  g = g ^ rot (c, 10);
  h = h ^ rot (a, 11);
  i = i ^ rot (b, 11);
  
  j = rot (*x, 21);
  k = rot (*y, 21);
  l = rot (*z, 21);
  
  j = j ^ a;
  k = k ^ b;
  l = l ^ c;
  
  g = g ^ rot (e, 26);
  h = h ^ rot (f, 26);
  i = i ^ rot (d, 27);
  
  a = rot (*y, 1) ^ rot (*z, 22);
  b = rot (*z, 1) ^ rot (*x, 23);
  c = rot (*x, 2) ^ rot (*y, 23);
  
  j = j ^ rot (c, 18);
  k = k ^ rot (a, 19);
  l = l ^ rot (b, 19);
  
  d = d ^ a ^ rot (d, 8);
  e = e ^ b ^ rot (e, 8);
  f = f ^ c ^ rot (f, 8);
  
  a = d ^ g ^ rot (e, 11);
  b = e ^ h ^ rot (f, 11);
  c = f ^ i ^ rot (d, 12);
  
  a = a ^ j ^ rot (a, 30);
  b = b ^ k ^ rot (b, 30);
  c = c ^ l ^ rot (c, 30);
  
  *x = a;
  *y = b;
  *z = c;
}

void print_interleaved(word a, word b, word c) {
  int tmp = 0;
  int tmpbits = 0;
  for (int i=31; i>=0; i--) {
    tmp = (tmp<<3) + 4*((a>>i)&1) + 2*((b>>i)&1) + 1*((c>>i)&1);
    tmpbits += 3;
    if (tmpbits >= 4) {
      int top4 = tmp>>(tmpbits-4);
      printf ("%01x", top4);
      tmp -= top4<<(tmpbits-4);
      tmpbits -= 4;
    }
  }
}

int main ()
{
  word x = 0;
  word y = 0;
  word z = 1;
  
  printf ("x = %u \ny = %u \nz = %u\n\n", x, y, z);
  
  M (&x, &y, &z);
  printf ("After the application of L_32x3 LBox:\n");
  printf ("x = %u \ny = %u \nz = %u\n\n", x, y, z);
  printf ("interleaved representation: ");
  print_interleaved(x, y, z);
  printf("\n\n");
  
  I (&x, &y, &z);
  printf ("After the application of L_32x3 LBox inverse:\n");
  printf ("x = %u \ny = %u \nz = %u\n", x, y, z);
  
  return 0;
}
