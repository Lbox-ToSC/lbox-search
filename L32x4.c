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


// Implementation of L_32x4 LBox
void M (word *x, word *y, word *z, word *t)
{
  word a, b, c, d, e, f, g, h, i, j, k, l;
  
  a = *x ^ rot (*y, 19) ^ rot (*t, 24);
  b = *y ^ rot (*z, 19) ^ rot (*x, 25);
  c = *z ^ rot (*t, 19) ^ rot (*y, 25);
  d = *t ^ rot (*x, 20) ^ rot (*z, 25);
  
  e = a ^ rot (c, 2);
  f = b ^ rot (d, 2);
  g = c ^ rot (a, 3);
  h = d ^ rot (b, 3);
  
  i = e ^ rot (h, 8);
  j = f ^ rot (e, 9);
  k = g ^ rot (f, 9);
  l = h ^ rot (g, 9);
  
  e = i ^ rot (d, 30);  
  f = j ^ rot (a, 31);
  g = k ^ rot (b, 31);
  h = l ^ rot (c, 31);
  
  a = e ^ rot (k, 3);
  b = f ^ rot (l, 3);
  c = g ^ rot (i, 4);
  d = h ^ rot (j, 4);
  
  *x = a;
  *y = b;
  *z = c;
  *t = d;
}

// Implementation of L_32x4 LBox inverse
void I (word *x, word *y, word *z, word *t)
{
  word a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, u;

  a = rot (*t, 2) ^ rot (*x, 28);
  b = rot (*x, 3) ^ rot (*y, 28);
  c = rot (*y, 3) ^ rot (*z, 28);
  d = rot (*z, 3) ^ rot (*t, 28);
  
  e = rot (*y, 15) ^ a;
  f = rot (*z, 15) ^ b;
  g = rot (*t, 15) ^ c;
  h = rot (*x, 16) ^ d;
  
  a = rot (*y, 25) ^ rot (a, 21);
  b = rot (*z, 25) ^ rot (b, 21);
  c = rot (*t, 25) ^ rot (c, 21);
  d = rot (*x, 26) ^ rot (d, 21);
  
  i = rot (*x, 9) ^ rot (*t, 19);
  j = rot (*y, 9) ^ rot (*x, 20);
  k = rot (*z, 9) ^ rot (*y, 20);
  l = rot (*t, 9) ^ rot (*z, 20);
  
  m = i ^ rot (k, 2) ^ rot (e, 7);
  n = j ^ rot (l, 2) ^ rot (f, 7);
  o = k ^ rot (i, 3) ^ rot (g, 7);
  p = l ^ rot (j, 3) ^ rot (h, 7);
  
  q = rot (*t, 3) ^ rot (*y, 7);
  r = rot (*x, 4) ^ rot (*z, 7);
  s = rot (*y, 4) ^ rot (*t, 7);
  u = rot (*z, 4) ^ rot (*x, 8);
  
  e = e ^ rot (q, 12);
  f = f ^ rot (r, 12);
  g = g ^ rot (s, 12);
  h = h ^ rot (u, 12);
  
  a = a ^ q ^ rot (r, 3);
  b = b ^ r ^ rot (s, 3);
  c = c ^ s ^ rot (u, 3);
  d = d ^ u ^ rot (q, 4);
  
  q = a ^ e ^ rot (d, 0);
  r = b ^ f ^ rot (a, 1);
  s = c ^ g ^ rot (b, 1);
  u = d ^ h ^ rot (c, 1);
  
  k = rot (*x, 1) ^ rot (*t, 26) ^ rot (k, 4);
  l = rot (*y, 1) ^ rot (*x, 27) ^ rot (l, 4);
  i = rot (*z, 1) ^ rot (*y, 27) ^ rot (i, 5);
  j = rot (*t, 1) ^ rot (*z, 27) ^ rot (j, 5);
  
  a = k ^ m ^ rot (m, 12);
  b = l ^ n ^ rot (n, 12);
  c = i ^ o ^ rot (o, 12);
  d = j ^ p ^ rot (p, 12);
  
  m = q ^ a ^ rot (c, 18);
  n = r ^ b ^ rot (d, 18);
  o = s ^ c ^ rot (a, 19);
  p = u ^ d ^ rot (b, 19);
  
  *x = m;
  *y = n;
  *z = o;
  *t = p;
}

void print_interleaved(word a, word b, word c, word d) {
  for (int i=31; i>=0; i--) {
    printf ("%01x", 8*((a>>i)&1) + 4*((b>>i)&1) + 2*((c>>i)&1) + 1*((d>>i)&1));
  }
}

int main ()
{
  word x = 0;
  word y = 0;
  word z = 0;
  word t = 1;

  printf ("x = %08x \ny = %08x \nz = %08x \nt = %08x\n\n", x, y, z, t);
  
  M (&x, &y, &z, &t);
  printf ("After the application of L_32x4 LBox:\n");
  printf ("x = %08x \ny = %08x \nz = %08x \nt = %08x\n", x, y, z, t);
  printf ("interleaved representation: ");
  print_interleaved(x, y, z, t);
  printf("\n\n");
  
  I (&x, &y, &z, &t);
  printf ("After the application of L_32x4 LBox inverse:\n");
  printf ("x = %08x \ny = %08x \nz = %08x \nt = %08x\n", x, y, z, t);
  
  return 0;
}
