/*
 * vf3.h
 *
 *  Created on: Jun 21, 2021
 *       Author: Gustavo Banegas
 */

#ifndef VF3_H_
#define VF3_H_

#define popcount __builtin_popcountll

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <immintrin.h>

#include "libpopcnt.h"

#include "definitions.h"

#include "f3.h"

#include "prng.h"

f3_vector f3_vector_new(size_t length);

f3_vector f3_vector_new_copy_from_array(const size_t length, const uint8_t *a);

size_t f3_vector_weight(const f3_vector *vec);

void f3_vector_free(f3_vector *x);

void f3_vector_zero(f3_vector *x);

void f3_vector_setcoeff(f3_vector *x, size_t j, uint8_t a);

uint8_t f3_vector_coeff_iszero(f3_vector *x, int j);

uint8_t f3_vector_coeff_isone(f3_vector *x, int j);

uint8_t f3_vector_coeff_istwo(f3_vector *x, int j);

uint8_t f3_vector_get_coeff(f3_vector *x, int j);

uint8_t f3_vector_mul_sum(const f3_vector *x, const f3_vector *y);

void f3_vector_print(const f3_vector *a);

void f3_vector_print_big(const f3_vector *a);

uint8_t f3_vector_sum(f3_vector *x);

uint8_t f3_vector_dotproduct(const f3_vector *x, const f3_vector *y);
uint8_t f3_vector_dotproduct_slow(const f3_vector *x, const f3_vector *y);

void f3_vector_mul(const f3_vector *x, const f3_vector *y, f3_vector *res);

void f3_vector_rand(f3_vector *a, prng_t *PRNG);

void f3_vector_set_from_array(f3_vector *x, const uint8_t *v,
		const size_t length);

uint8_t f3_vector_cmp(const f3_vector *a, const f3_vector *b, size_t nr_bits);

static inline void f3_vector_sum_inplace(f3_vector *a, const f3_vector *b) {
	for (size_t i = 0; i < a->alloc; i++) {
		//f3_inplace_sum(&a->x[i], &b->x[i]);

#ifdef SIXOP
		wave_word tmp = 0;
		tmp = (~(a->r0[i] ^ b->r0[i])) & ((a->r1[i] ^ b->r1[i]) ^ a->r0[i]);
		a->r0[i] = (~(a->r1[i] ^ b->r1[i]))
				& ((a->r0[i] ^ b->r0[i]) ^ a->r1[i]);
		a->r1[i] = tmp;

#elif SEVENOP

 #endif
		/*	vf3_add_aux(&a->x[i].r0, &a->x[i].r1, a->x[i].r0, a->x[i].r1,
		 b->x[i].r0, b->x[i].r1);*/
	}
}

static inline void f3_vector_add(const f3_vector *a, const f3_vector *b,
		f3_vector *c) {
	for (size_t i = 0; i < a->alloc; i++) {
		//f3_inplace_sum(&a->x[i], &b->x[i]);

#ifdef SIXOP
		c->r1[i] = (~(a->r0[i] ^ b->r0[i]))
				& ((a->r1[i] ^ b->r1[i]) ^ a->r0[i]);
		c->r0[i] = (~(a->r1[i] ^ b->r1[i]))
				& ((a->r0[i] ^ b->r0[i]) ^ a->r1[i]);

#elif SEVENOP

 #endif
		/*	vf3_add_aux(&a->x[i].r0, &a->x[i].r1, a->x[i].r0, a->x[i].r1,
		 b->x[i].r0, b->x[i].r1);*/
	}
}
static inline void f3_vector_neg(const f3_vector *a, f3_vector *b) {
	for (int i = 0; i < a->alloc; i++) {
		b->r0[i] = a->r1[i];
		b->r1[i] = a->r0[i];
	}

}

// c = a - b
static inline void f3_vector_sub(const f3_vector *a, const f3_vector *b,
		f3_vector *c) {
	for (size_t i = 0; i < a->alloc; i++) {
		c->r0[i] = ((a->r0[i] ^ b->r1[i]) & (a->r1[i] ^ b->r0[i]))
				^ (a->r1[i] & b->r0[i]) ^ (a->r0[i] ^ b->r1[i]);
		c->r1[i] = ((a->r0[i] ^ b->r1[i]) & (a->r1[i] ^ b->r0[i]))
				^ (a->r1[i] ^ b->r0[i]) ^ (a->r0[i] & b->r1[i]);

	}
}

static inline void f3_vector_sub_inplace(f3_vector *a, const f3_vector *b) {
	for (size_t i = 0; i < a->alloc; i++) {
		wave_word tmp = ((a->r0[i] ^ b->r1[i]) & (a->r1[i] ^ b->r0[i]))
				^ (a->r1[i] & b->r0[i]) ^ (a->r0[i] ^ b->r1[i]);
		a->r1[i] = ((a->r0[i] ^ b->r1[i]) & (a->r1[i] ^ b->r0[i]))
				^ (a->r1[i] ^ b->r0[i]) ^ (a->r0[i] & b->r1[i]);
		a->r0[i] = tmp;

	}
}

static inline uint8_t f3_vector_coeff(const f3_vector *x, int j) {
	int i;
	i = j / WORD_LENGTH;
	j %= WORD_LENGTH;
	uint8_t c = ((((x->r0[i] >> j) & 1) - ((x->r1[i] >> j) & 1)) + 3) % 3;
	return c;
}

f3_vector f3_vector_one_new(size_t length);

static inline void f3_vector_add_avx(const f3_vector *a, const f3_vector *b,
		f3_vector *c) {

	for (int offset = 0; offset < a->alloc; offset = offset + 4) {
		double *addr_a = (double*) (&a->r0[offset]);
		__m256d a_l = _mm256_loadu_pd(addr_a);
		__m256d a_h = _mm256_loadu_pd((double*) a->r1 + offset);
		__m256d b_l = _mm256_loadu_pd((double*) b->r0 + offset);
		__m256d b_h = _mm256_loadu_pd((double*) b->r1 + offset);

		__m256d res_1 = _mm256_xor_pd(a_l, b_l);
		__m256d res_2 = _mm256_xor_pd(a_h, b_h);
		__m256d res_3 = _mm256_xor_pd(res_2, a_l);
		__m256d c_high = _mm256_andnot_pd(res_1, res_3);
		_mm256_storeu_pd((double*) c->r1 + offset, c_high);

		res_1 = _mm256_xor_pd(a_h, b_h);
		res_2 = _mm256_xor_pd(a_l, b_l);
		res_3 = _mm256_xor_pd(res_2, a_h);
		__m256d c_low = _mm256_andnot_pd(res_1, res_3);
		//__m256i* ptr_l = (__m256i*) (c->r0 + ptr);
		_mm256_storeu_pd((double*) c->r0 + offset, c_low);

	}

	//__m256i res = _mm256_xor_epi64(a->r0, b->r0);
	/*
	 for (size_t i = 0; i < a->alloc; i++) {
	 //f3_inplace_sum(&a->x[i], &b->x[i]);

	 #ifdef SIXOP
	 c->r1[i] = (~(a->r0[i] ^ b->r0[i]))
	 & ((a->r1[i] ^ b->r1[i]) ^ a->r0[i]);

	 c->r0[i] = (~(a->r1[i] ^ b->r1[i]))
	 & ((a->r0[i] ^ b->r0[i]) ^ a->r1[i]);

	 #elif SEVENOP

	 #endif
	 vf3_add_aux(&a->x[i].r0, &a->x[i].r1, a->x[i].r0, a->x[i].r1, b->x[i].r0,
	 b->x[i].r1);
	 }*/
}

static inline void f3_vector_sum_inplace_avx(const f3_vector *a,
		const f3_vector *b) {

	int offset = 0;
	for (offset = 0; offset < a->alloc - 4; offset = offset + 4) {
		double *addr_a = (double*) (&a->r0[offset]);
		__m256d a_l = _mm256_loadu_pd(addr_a);
		__m256d a_h = _mm256_loadu_pd((double*) a->r1 + offset);
		__m256d b_l = _mm256_loadu_pd((double*) b->r0 + offset);
		__m256d b_h = _mm256_loadu_pd((double*) b->r1 + offset);

		__m256d res_1 = _mm256_xor_pd(a_l, b_l);
		__m256d res_2 = _mm256_xor_pd(a_h, b_h);
		__m256d res_3 = _mm256_xor_pd(res_2, a_l);
		__m256d c_high = _mm256_andnot_pd(res_1, res_3);
		_mm256_storeu_pd((double*) a->r1 + offset, c_high);

		res_1 = _mm256_xor_pd(a_h, b_h);
		res_2 = _mm256_xor_pd(a_l, b_l);
		res_3 = _mm256_xor_pd(res_2, a_h);
		__m256d c_low = _mm256_andnot_pd(res_1, res_3);
		//__m256i* ptr_l = (__m256i*) (c->r0 + ptr);
		_mm256_storeu_pd((double*) a->r0 + offset, c_low);

	}
	if (offset < a->alloc) {
		if (a->alloc % offset == 2) {
			double *addr_a = (double*) (&a->r0[offset]);
			__m128d a_l = _mm_load_pd(addr_a);
			__m128d a_h = _mm_load_pd((double*) a->r1 + offset);
			__m128d b_l = _mm_load_pd((double*) b->r0 + offset);
			__m128d b_h = _mm_load_pd((double*) b->r1 + offset);

			__m128d res_1 = _mm_xor_pd(a_l, b_l);
			__m128d res_2 = _mm_xor_pd(a_h, b_h);
			__m128d res_3 = _mm_xor_pd(res_2, a_l);
			__m128d c_high = _mm_andnot_pd(res_1, res_3);
			_mm_storeu_pd((double*) a->r1 + offset, c_high);

			res_1 = _mm_xor_pd(a_h, b_h);
			res_2 = _mm_xor_pd(a_l, b_l);
			res_3 = _mm_xor_pd(res_2, a_h);
			__m128d c_low = _mm_andnot_pd(res_1, res_3);
			//__m256i* ptr_l = (__m256i*) (c->r0 + ptr);
			_mm_storeu_pd((double*) a->r0 + offset, c_low);
		} else if (a->alloc % offset == 1) {
			wave_word tmp = 0;
			tmp = (~(a->r0[offset] ^ b->r0[offset]))
					& ((a->r1[offset] ^ b->r1[offset]) ^ a->r0[offset]);
			a->r0[offset] = (~(a->r1[offset] ^ b->r1[offset]))
					& ((a->r0[offset] ^ b->r0[offset]) ^ a->r1[offset]);
			a->r1[offset] = tmp;
		}
	}
}

static inline void f3_vector_sub_avx(const f3_vector *a, const f3_vector *b,
		f3_vector *c) {
	int offset = 0;
	for (offset = 0; offset < a->alloc - 4; offset = offset + 4) {
		__m256d a_l = _mm256_loadu_pd((double*) a->r0 + offset);
		__m256d a_h = _mm256_loadu_pd((double*) a->r1 + offset);

		__m256d b_l = _mm256_loadu_pd((double*) b->r0 + offset);
		__m256d b_h = _mm256_loadu_pd((double*) b->r1 + offset);

		__m256d res_1 = _mm256_xor_pd(a_l, b_h);

		__m256d res_2 = _mm256_xor_pd(a_h, b_l);

		__m256d res_3 = _mm256_and_pd(res_1, res_2);

		__m256d res_4 = _mm256_and_pd(a_h, b_l);

		__m256d res_5 = _mm256_xor_pd(res_3, res_4);

		__m256d c_low = _mm256_xor_pd(res_5, res_1);

		_mm256_storeu_pd((double*) c->r0 + offset, c_low);

		res_4 = _mm256_and_pd(a_l, b_h);
		res_5 = _mm256_xor_pd(a_h, b_l);
		c_low = _mm256_xor_pd(res_4, res_5);
		__m256d c_high = _mm256_xor_pd(res_3, c_low);

		_mm256_storeu_pd((double*) c->r1 + offset, c_high);

	}

	if (offset < a->alloc) {
		__m128d a_l = _mm_load_pd((double*) a->r0 + offset);
		__m128d a_h = _mm_load_pd((double*) a->r1 + offset);
		__m128d b_l = _mm_load_pd((double*) b->r0 + offset);
		__m128d b_h = _mm_load_pd((double*) b->r1 + offset);

		__m128d res_1 = _mm_xor_pd(a_l, b_h);

		__m128d res_2 = _mm_xor_pd(a_h, b_l);

		__m128d res_3 = _mm_and_pd(res_1, res_2);

		__m128d res_4 = _mm_and_pd(a_h, b_l);

		__m128d res_5 = _mm_xor_pd(res_3, res_4);

		__m128d c_low = _mm_xor_pd(res_5, res_1);

		_mm_storeu_pd((double*) c->r0 + offset, c_low);

		res_4 = _mm_and_pd(a_l, b_h);
		res_5 = _mm_xor_pd(a_h, b_l);
		c_low = _mm_xor_pd(res_4, res_5);
		__m128d c_high = _mm_xor_pd(res_3, c_low);
		//__m256i* ptr_l = (__m256i*) (c->r0 + ptr);
		_mm_storeu_pd((double*) c->r1 + offset, c_high);
	}
}

static inline void f3_vector_sub_inplace_avx(f3_vector *a, const f3_vector *b) {
	int offset = 0;
	for (offset = 0; offset < a->alloc - 4; offset = offset + 4) {
		__m256d a_l = _mm256_loadu_pd((double*) a->r0 + offset);
		__m256d a_h = _mm256_loadu_pd((double*) a->r1 + offset);

		__m256d b_l = _mm256_loadu_pd((double*) b->r0 + offset);
		__m256d b_h = _mm256_loadu_pd((double*) b->r1 + offset);

		__m256d res_1 = _mm256_xor_pd(a_l, b_h);

		__m256d res_2 = _mm256_xor_pd(a_h, b_l);

		__m256d res_3 = _mm256_and_pd(res_1, res_2);

		__m256d res_4 = _mm256_and_pd(a_h, b_l);

		__m256d res_5 = _mm256_xor_pd(res_3, res_4);

		__m256d c_low = _mm256_xor_pd(res_5, res_1);

		_mm256_storeu_pd((double*) a->r0 + offset, c_low);

		res_4 = _mm256_and_pd(a_l, b_h);
		res_5 = _mm256_xor_pd(a_h, b_l);
		c_low = _mm256_xor_pd(res_4, res_5);
		__m256d c_high = _mm256_xor_pd(res_3, c_low);

		_mm256_storeu_pd((double*) a->r1 + offset, c_high);

	}

	if (offset < a->alloc) {

		if (a->alloc - offset == 2) {
			__m128d a_l = _mm_load_pd((double*) a->r0 + offset);
			__m128d a_h = _mm_load_pd((double*) a->r1 + offset);
			__m128d b_l = _mm_load_pd((double*) b->r0 + offset);
			__m128d b_h = _mm_load_pd((double*) b->r1 + offset);

			__m128d res_1 = _mm_xor_pd(a_l, b_h);

			__m128d res_2 = _mm_xor_pd(a_h, b_l);

			__m128d res_3 = _mm_and_pd(res_1, res_2);

			__m128d res_4 = _mm_and_pd(a_h, b_l);

			__m128d res_5 = _mm_xor_pd(res_3, res_4);

			__m128d c_low = _mm_xor_pd(res_5, res_1);

			_mm_storeu_pd((double*) a->r0 + offset, c_low);

			res_4 = _mm_and_pd(a_l, b_h);
			res_5 = _mm_xor_pd(a_h, b_l);
			c_low = _mm_xor_pd(res_4, res_5);
			__m128d c_high = _mm_xor_pd(res_3, c_low);
			//__m256i* ptr_l = (__m256i*) (c->r0 + ptr);
			_mm_storeu_pd((double*) a->r1 + offset, c_high);
		}
	}
}

#endif /* VF3_H_ */
