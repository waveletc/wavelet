/*
 ============================================================================
 Name        : wave.c
 Author      : Gustavo Banegas
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <inttypes.h>

#include "sage.h"

#include "definitions.h"
#include "params.h"
#include "commons.h"

#include "f3.h"
#include "vf3.h"
#include "mf3.h"

#include "keygen.h"
#include "sign.h"
#include "verify.h"
#include "compress.h"
#include "file_save.h"

#include "commons.h"

#include "cpucycles.h"

#define MAXMLEN 16

#define ITERATIONS 10

void swap(unsigned long long *a, unsigned long long *b) {
	unsigned long long t = *a;
	*a = *b;
	*b = t;
}

// A function to implement bubble sort
void bubbleSort(unsigned long long arr[], int n) {
	int i, j;
	for (i = 0; i < n - 1; i++)

		// Last i elements are already in place
		for (j = 0; j < n - i - 1; j++)
			if (arr[j] > arr[j + 1])
				swap(&arr[j], &arr[j + 1]);
}

void print_message(uint8_t *m, size_t len) {
	printf("{ ");
	for (size_t i = 0; i < len; i++) {
		printf("%.2x, ", m[i]);
	}
	printf("} \n");
}

void print_signature(f3_vector *sig) {

	for (int i = 0; i < sig->alloc; i++) {
		printf("0x%" PRIx64 ", ", sig->r0[i]);
	}
	printf("\n");

	for (int i = 0; i < sig->alloc; i++) {
		printf("0x%" PRIx64 ", ", sig->r1[i]);
	}
	printf("\n");
}

unsigned long long time_sign_f[ITERATIONS];
unsigned long long time_verify_f[ITERATIONS];
unsigned long long time_stverify_f[ITERATIONS];
unsigned long long time_compress_f[ITERATIONS];
unsigned long long time_decompress_f[ITERATIONS];

void generate_t(mf3 *t) {
	int sign = 0;
	for (int i = 0; i < t->n_rows - 1; i++) {
		for (int j = 0; j < t->n_cols - 1; j++) {
			if (i == j) {
				if (sign) {
					mf3_setcoeff(t, i, j, 2);
					sign = 0;
				} else {
					mf3_setcoeff(t, i, j, 1);
					mf3_setcoeff(t, i + 1, j, 1);
					mf3_setcoeff(t, i, j + 1, 1);
					sign = 1;
				}
			}
		}
	}

}

int main(void) {

	init(PARAMS_ID);
	unsigned long long start, end;
	memset(time_sign_f, 0, sizeof(unsigned long long) * ITERATIONS);
	memset(time_verify_f, 0, sizeof(unsigned long long) * ITERATIONS);
	memset(time_stverify_f, 0, sizeof(unsigned long long) * ITERATIONS);
	memset(time_compress_f, 0, sizeof(unsigned long long) * ITERATIONS);
	memset(time_decompress_f, 0, sizeof(unsigned long long) * ITERATIONS);

	wave_pk_t pk;
	wave_sk_t sk;

	f3_vector signature = f3_vector_new(N);

	f3_vector m_hash = f3_vector_new(N - K);
	uint8_t mi[MAXMLEN] = { 0 };

	printf("keygen starting....\n");
	start = cpucycles();
	keygen(&sk, &pk);
	end = cpucycles();
	unsigned long long cpu_time_used = (end - start);

	for (int i = 0; i < ITERATIONS; i++) {
		printf("Iteration: %d\n", i);
		f3_vector_zero(&signature);
		f3_vector_zero(&m_hash);
		uint8_t sig_compressed[N] = { 0 };
		size_t size_out = 0;
		memset(mi, 0, MAXMLEN * sizeof(uint8_t));
		randombytes(mi, MAXMLEN);

		uint8_t salt[SALT_SIZE] = { 0 };
		start = cpucycles();
		sign(&signature, &m_hash, salt, mi, MAXMLEN, sk);
		end = cpucycles();
		time_sign_f[i] = (end - start);

		start = cpucycles();
		int res = verify(salt, mi, MAXMLEN, &signature, pk);
		end = cpucycles();
		time_verify_f[i] = (end - start);
		if (res != 1) {
			printf("verification failed for message seed %d!\n", 0);
			return -1;
		}

		res = 0;

		f3_vector a = f3_vector_new(K);
		memcpy(a.r0, signature.r0 + 45, 88 * sizeof(wave_word));
		memcpy(a.r1, signature.r1 + 45, 88 * sizeof(wave_word));
		shift7(a.r0, a.alloc);
		shift7(a.r1, a.alloc);

		//print_signature(&a);
		start = cpucycles();
		encode_sig_triple(&a, a.size, sig_compressed, &size_out);
		end = cpucycles();
		time_compress_f[i] = (end - start);

		//printf("size: %zu\n", size_out);
		/*for (int j = 0; j < size_out; j++) {
		 printf("%x", sig_compressed[j]);
		 }
		 printf("\n");*/

		//cpu_time_used_comp += (end - start);
		f3_vector signature_out = f3_vector_new(K);
		start = cpucycles();
		decode_sig_triple(sig_compressed, size_out, &signature_out);
		end = cpucycles();
		time_decompress_f[i] = (end - start);

		int cmp_t = memcmp(signature_out.r0, a.r0, a.alloc * sizeof(wave_word));
		if (cmp_t != 0) {
			printf("fail: r0\n");
			return -1;
		}
		cmp_t = memcmp(signature_out.r1, a.r1, a.alloc * sizeof(wave_word));
		if (cmp_t != 0) {
			printf("fail: r1\n");
			return -1;
		}

		f3_vector one = f3_vector_one_new(N);
		f3_vector u = f3_vector_new(N - K);
		mf3 *r_t = mf3_transpose(pk);
		mf3_mv_mul_v(pk, &one, &u);

		start = cpucycles();
		//c_start = cpucycles();
		res = STverify(salt, mi, MAXMLEN, &a, sig_compressed, size_out, &u,
				r_t);
		//c_final = cpucycles();
		end = cpucycles();
		time_stverify_f[i] = (end - start);
		if (res != 1) {
			printf("verification failed for message seed %d!\n", 0);
			return -1;
		}

		mf3 *m_t = mf3_new(r_t->n_rows, r_t->n_cols);
		for (int a = 0; a < (K - 1) / 2; a++) {
			f3_vector_add(&r_t->row[2 * a], &r_t->row[(2 * a) + 1],
					&m_t->row[2 * a]);
			f3_vector_sub(&r_t->row[2 * a], &r_t->row[2 * a + 1],
					&m_t->row[2 * a + 1]);

		}
		f3_vector_neg(&r_t->row[K - 1], &m_t->row[K - 1]);

		big_save_public_key(m_t);
		big_save_m_hash_vector(&m_hash, i);
		print_signature(&a);
		//big_save_truncated_signature_vector(&a, i);

		start = cpucycles();
		//c_start = cpucycles();
		res = Nverify(salt, mi, MAXMLEN, &a, m_t);
		//c_final = cpucycles();
		end = cpucycles();
		if (res != 1) {
			printf("Nverify failed for message seed %d!\n", 0);
			return -1;
		}

		f3_vector_free(&one);
		f3_vector_free(&a);
		f3_vector_free(&u);
		f3_vector_free(&signature_out);
		mf3_free(r_t);
		mf3_free(m_t);

	}
	bubbleSort(time_stverify_f, ITERATIONS - 1);
	bubbleSort(time_verify_f, ITERATIONS - 1);
	bubbleSort(time_sign_f, ITERATIONS - 1);
	bubbleSort(time_decompress_f, ITERATIONS - 1);
	bubbleSort(time_compress_f, ITERATIONS - 1);

	printf("keygen() took %llu cycles to execute \n", cpu_time_used);
	printf("sign() took %llu cycles to execute \n",
			time_sign_f[ITERATIONS / 2]);
	printf("verify() took %llu cycles to execute \n",
			time_verify_f[ITERATIONS / 2]);
	printf("STverify() took %llu cycles to execute \n",
			time_stverify_f[ITERATIONS / 2]);

	printf("compression: took %llu cycles to execute\n",
			time_compress_f[ITERATIONS / 2]);

	printf("decompression: took %llu cycles to execute\n",
			time_decompress_f[ITERATIONS / 2]);

	wave_sk_clear(&sk);
	printf("cleared SK...\n");
	wave_pk_clear(pk);
	printf("cleared PK...\n");
	f3_vector_free(&signature);
	f3_vector_free(&m_hash);

	cleanup();

	return EXIT_SUCCESS;
}
