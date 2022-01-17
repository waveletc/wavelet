/*
 * f3.h
 *
 *  Created on: Jun 21, 2021
 *       Author: Gustavo Banegas
 */

#ifndef F3_H_
#define F3_H_

#include <stdio.h>
#include <stdlib.h>

#include "prng.h"

#include "definitions.h"



#define f3_mul(a, b) (((a) * (b))%3)
#define f3_add(a, b) (((a) + (b))%3)

#define f3_neg(a) (((-(a))+3)%3)

/*
 * a <- a + b
 */
static inline void f3_inplace_sum(f3 *a, const f3 *b) {

#ifdef SIXOP
	wave_word tmp = 0;
	tmp = (~(a->r0 ^ b->r0)) & ((a->r1 ^ b->r1) ^ a->r0);
	a->r0 = (~(a->r1 ^ b->r1)) & ((a->r0 ^ b->r0) ^ a->r1);
	a->r1 = tmp;

#elif SEVENOP

#endif
}
;

/*
 * c <- a + b
 */
static inline void f3_sum(const f3 *a, const f3 *b, f3 *c) {

#ifdef SIXOP
	c->r1 = (~(a->r0 ^ b->r0)) & ((a->r1 ^ b->r1) ^ a->r0);
	c->r0 = (~(a->r1 ^ b->r1)) & ((a->r0 ^ b->r0) ^ a->r1);
#elif SEVENOP

#endif

}
;



void f3_print(const f3 *a);
void f3_compact_print(const wave_word r0, const wave_word r1);

wave_word f3_array_weight(uint8_t *v, int length);

uint8_t f3_rand(prng_t *PRNG);

uint8_t f3_randnonzero(prng_t *PRNG);

#endif /* F3_H_ */
