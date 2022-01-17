/*
 * f3.c
 *
 *  Created on: Jun 21, 2021
 *       Author: Gustavo Banegas
 */

#include "f3.h"
#include "libpopcnt.h"



void f3_print(const f3 *a) {
	wave_word one = 1;
	for (size_t i = 0; i < WORD_LENGTH; i++) {

		wave_word pos = one << i;
#ifdef SIXOP
		if (a->r0 & pos) {
			printf("1, ");
		} else if (a->r1 & pos) {
			printf("2, ");
		} else {
			printf("0, ");
		}
#elif SEVENOP

#endif
	}
	printf("\n");

}





void f3_compact_print(const wave_word r0, const wave_word r1) {
	wave_word one = 1;
	for (size_t i = 0; i < WORD_LENGTH; i++) {
		wave_word pos = one << i;
#ifdef SIXOP
		if (r0 & pos) {
			printf("1,");
		} else if (r1 & pos) {
			printf("2,");
		} else {
			printf("0,");
		}
#elif SEVENOP

	#endif
	}
	//printf("\n");
}

inline uint8_t f3_randnonzero(prng_t *PRNG) {
	int8_t c = 1 - (2 * rnd_bit(PRNG)); // 0 or 1
	c = (c + 3) % 3;
	return c; // 1 or -1
}

uint8_t f3_rand(prng_t *PRNG) {
	uint8_t c = rnd_trit(PRNG); // 0, 1, or 2
	return c; // 1, 0, or -1
}

// Hamming weight of v
wave_word f3_array_weight(uint8_t *v, int length) {
	//int i;
	wave_word j = 0;
	j = popcnt(v, length*sizeof(uint8_t));
	return j;
}
