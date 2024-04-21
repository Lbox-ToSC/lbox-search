# Design of a Linear Layer Optimised for Bitsliced 32-bit Implementation

This repository contains additional data and code from the paper:

> [*Design of a Linear Layer Optimised for Bitsliced 32-bit Implementation*](https://tosc.iacr.org/index.php/ToSC/article/view/11412)
> Gaëtan Leurent and Clara Pernot
> ToSC 2024(1)



## Repository content

This repository includes:

1. Code to search for efficient linear transforms, in file
   `find_lbox.c`;
   
2. Code to search for a relatively efficient implementation of a linear
   transform, in file `heuristic_implem.c`;

3. Reference implementation of the two linear transforms described in
   the paper, in files `L32x3.c` and `L32x4.c`.

## 1. Finding efficient linear transforms: `find_lbox.c`

This program implements the main search of the paper.
There are three ways t use the program

1. Look for random circulant matrices with high branch number
```
./find_lbox -r
```

2. Look for circulant matrices with efficient implementation and high branch number
```
./find_lbox
```

3. Verify the branch number of a given LBox (the example corresponds to L32x4 from the paper)
```
./find_lbox 44032028488021000802200841114541
```


Important parameters are `#define`ed at the top of the file:

- `WORD_SIZE`: parameter $\ell$ in the paper
- `WORD_NB`: parameter $w$ in the paper
- `TARGET_WEIGHT`: maximum weight for the search (inclusive).  In order to proove that the branch number is at least w, set `TARGET_WEIGHT` to w-1 in order to verify that no words of weight w or less are found.
- `ITER`: number of iterations (log2)
- `MAX_XOR`: number of rotation+XOR steps in the implementation
- `COMBINE3ROWS`: use combinaisons of up to three rows ($d=3$ in the paper). Otherwise, use combinaisons of up to two rows ($d=2$)

Some parameter examples follow

### Parameters for an LBox over 3 32-bit words

```
#define WORD_SIZE 32
#define WORD_NB 3
#define TARGET_WEIGHT 18
#define ITER 10
#define MAX_XOR 6
#define COMBINE3ROWS
```

### Parameters for an LBox over 4 32-bit words

```
#define WORD_SIZE 32
#define WORD_NB 4
#define TARGET_WEIGHT 20
#define ITER 10
#define MAX_XOR 6
```

### Parameters for an LBox over 4 16-bit words

```
#define WORD_SIZE 16
#define WORD_NB 4
#define TARGET_WEIGHT 11
#define ITER 10
#define MAX_XOR 5
```

### Parameters for an LBox over 8 16-bit words

```
#define WORD_SIZE 16
#define WORD_NB 4
#define TARGET_WEIGHT 13
#define ITER 12
#define MAX_XOR 6
```


## 2. Finding an implementation of a linear transforms: `heuristic_implem.c`

This code uses the `_BitInt` feature of C23 for variable-length
integers.  Therefore you need a recent version of `clang` (>=15) in
order to compile it.

The word size is assumed to bit 32 ($\ell = 32$).
The `WIDTH` parameter defined at the top of the file corresponds to the parameter $w$ in the paper

Example (the example corresponds to the inverse of L32x4 from the paper):

```
./heuristic_implem 81034926c0d0e290b0606d93f2049851
```
