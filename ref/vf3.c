/*
 * vf3.c
 *
 *  Created on: Jun 21, 2021
 *       Author: Gustavo Banegas
 */
#include <popcntintrin.h>

#include <nmmintrin.h>

#include "vf3.h"

f3_vector f3_vector_new(size_t length) {

	f3_vector x;
	x.size = length;
	x.alloc = 1 + ((length - 1) / WORD_LENGTH);

	x.r0 = (wave_word*) calloc(x.alloc, sizeof(wave_word));
	x.r1 = (wave_word*) calloc(x.alloc, sizeof(wave_word));

	return x;
}

f3_vector f3_vector_new_copy_from_array(const size_t length, const uint8_t *v) {

	f3_vector x;
	x.size = length;
	x.alloc = 1 + ((length - 1) / WORD_LENGTH);

	x.r0 = (wave_word*) calloc(x.alloc, sizeof(wave_word));
	x.r1 = (wave_word*) calloc(x.alloc, sizeof(wave_word));

	size_t i, j;
	wave_word z;

	for (i = 0, j = 0, z = 1; j < length; j++, z <<= 1) {
		if (z == 0) {
			z = 1;
			i++;
		}
		if (v[j] == 1)
			x.r0[i] ^= z;
		else if (v[j] == 2)
			x.r1[i] ^= z;
	}

	return x;
}




void f3_vector_free(f3_vector *x) {
	free(x->r0);
	free(x->r1);
}

void f3_vector_zero(f3_vector *x) {
	memset(x->r0, 0, x->alloc * sizeof(wave_word));
	memset(x->r1, 0, x->alloc * sizeof(wave_word));
}

unsigned char compact(uint8_t *d, int size) {
	unsigned char c = 0;
	for (int i = 0; i < size; i++) {

		c *= 3;
		c += d[i];
	}

	return c;
}

void toarray(const f3_vector *x, uint8_t *a, int length) {
	int i, j;
	uint_fast64_t z;

	for (i = 0, j = 0, z = 1; j < length; j++, z <<= 1) {
		if (z == 0) {
			z = 1;
			i++;
		}
		if (x->r0[i] & z)
			a[j] = 1;
		else if (x->r1[i] & z)
			a[j] = 2;
		else
			a[j] = 0;
	}
}

void f3_vector_print(const f3_vector *a) {
	uint8_t *v = malloc(a->size);
	toarray(a, v, a->size);
	int i;
	for (i = 0; i < a->size - 5; i += 5) {
		printf("%x", compact(v + i, 5));
	}
	if (i < a->size) {

		printf("%x", compact(v + i, a->size - i));

	}
	printf("\n");
	free(v);
}

void f3_vector_print_big(const f3_vector *a) {
	for (size_t i = 0; i < a->alloc; i++) {
		f3_compact_print(a->r0[i], a->r1[i]);
		//printf(" | ");
	}
}

int binary_weight(const wave_word *x, size_t n) {
	size_t i, w = 0;
	for (i = 0; i < n; i++)
		w += popcount(x[i]);
	return w;
}

size_t f3_vector_weight(const f3_vector *vec) {
	size_t w = 0;
	w = popcnt(vec->r0, vec->alloc * 8) + popcnt(vec->r1, vec->alloc * 8);
	return w;
}

uint8_t f3_vector_sum(f3_vector *x) {
// int8_t = {0,1,-1}
	uint8_t j = (popcnt(x->r0, x->alloc * 8) + 2 * popcnt(x->r1, x->alloc * 8))
			% 3;
	return j;
}

f3_vector f3_vector_one_new(size_t length) {

	f3_vector x;
	x.size = length;
	x.alloc = 1 + ((length - 1) / WORD_LENGTH);

	x.r0 = (wave_word*) calloc(x.alloc, sizeof(wave_word));
	x.r1 = (wave_word*) calloc(x.alloc, sizeof(wave_word));

	memset(x.r0, 0xFFFFFFFF, x.alloc*sizeof(wave_word));

	return x;
}

/*f3_vector_mul(const f3_vector *x, const f3_vector *y, f3_vector *res) {
 for (size_t i = 0; i < x->alloc; i++) {

 res->r0[i] = (x->r0[i] & y->r0[i]) ^ (x->r1[i] & y->r1[i]);
 res->r1[i] = (x->r0[i] & y->r1[i]) ^ (x->r1[i] & y->r0[i]);

 }
 }*/
// z <- x * y
// componentwise product
//TODO: use AVX
uint8_t f3_vector_mul_sum(const f3_vector *x, const f3_vector *y) {
	size_t i;
	uint8_t j = 0;
	for (i = 0; i < x->alloc; i++) {
		wave_word r0, r1;
		r0 = (x->r0[i] & y->r0[i]) ^ (x->r1[i] & y->r1[i]);
		r1 = (x->r0[i] & y->r1[i]) ^ (x->r1[i] & y->r0[i]);
		j += ((popcnt(&r0, 8) + 2 * popcnt(&r1, 8)) % 3);
	}
	return (j % 3);
}


uint8_t f3_vector_dotproduct_slow(const f3_vector *x, const f3_vector *y) {
	f3_vector z = f3_vector_new(x->size);
	f3_vector_mul(x, y, &z);
	uint8_t a = f3_vector_sum(&z);
	f3_vector_free(&z);
	return a;
}

uint8_t f3_vector_dotproduct(const f3_vector *x, const f3_vector *y) {
//f3_vector z = f3_vector_new(x->size);
	uint8_t a = f3_vector_mul_sum(x, y);
//uint8_t a = f3_vector_sum(&z);
//f3_vector_free(&z);
	return a;
}

uint8_t f3_vector_get_coeff(f3_vector *x, int j) {
	if (f3_vector_coeff_iszero(x, j)) {
		return 0;
	} else if (f3_vector_coeff_isone(x, j)) {
		return 1;
	} else {
		return 2;
	}
}

uint8_t f3_vector_coeff_iszero(f3_vector *x, int j) {
	size_t i;
	i = j / WORD_LENGTH;
	j %= WORD_LENGTH;
	return (((x->r0[i] | x->r1[i]) >> j) & 1) == 0;
}

uint8_t f3_vector_coeff_isone(f3_vector *x, int j) {
	uint32_t i;
	i = j / WORD_LENGTH;
	j %= WORD_LENGTH;
	return (x->r0[i] >> j) & 1;
}

uint8_t f3_vector_coeff_istwo(f3_vector *x, int j) {
	uint32_t i;
	i = j / WORD_LENGTH;
	j %= WORD_LENGTH;
	return (x->r1[i] >> j) & 1;
}

void f3_vector_aux_multiplication(wave_word *r0, wave_word *r1, wave_word x0,
		wave_word x1, wave_word y0, wave_word y1) {
	wave_word a, b, c, d;
	a = x0 & y0;
	b = x1 & y1;
	c = x0 & y1;
	d = x1 & y0;
	*r0 = a ^ b;
	*r1 = c ^ d;
}

// z <- x * y
// componentwise product
void f3_vector_mul(const f3_vector *x, const f3_vector *y, f3_vector *res) {
	size_t i;
	for (i = 0; i < x->alloc; i++) {
		res->r0[i] = (x->r0[i] & y->r0[i]) ^ (x->r1[i] & y->r1[i]);
		res->r1[i] = (x->r0[i] & y->r1[i]) ^ (x->r1[i] & y->r0[i]);
	}
}

void f3_vector_setcoeff(f3_vector *x, size_t j, uint8_t a) {
	size_t i;
	i = j / WORD_LENGTH;
	j %= WORD_LENGTH;
	wave_word z = ((wave_word) 1) << j;
	if (a == 1) {
		x->r0[i] |= z;
		x->r1[i] &= ~z;
	} else if (a == 2) {
		x->r0[i] &= ~z;
		x->r1[i] |= z;
	} else {
		x->r0[i] &= ~z;
		x->r1[i] &= ~z;
	}
}

void f3_vector_rand(f3_vector *a, prng_t *PRNG) {

	for (size_t j = 0; j < a->size; j++) {
		f3_vector_setcoeff(a, j, rnd_trit(PRNG));
	}
}

void f3_vector_set_from_array(f3_vector *x, const uint8_t *v,
		const size_t length) {
	size_t i, j;
	wave_word z;

	f3_vector_zero(x);
	for (i = 0, j = 0, z = 1; j < length; j++, z <<= 1) {
		if (z == 0) {
			z = 1;
			i++;
		}
		if (v[j] == 1)
			x->r0[i] ^= z;
		else if (v[j] == 2)
			x->r1[i] ^= z;
	}
}
