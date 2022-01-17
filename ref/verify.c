#include "verify.h"
#include <time.h>

int verify(const uint8_t *salt, const uint8_t *message, const size_t mlen,
		const f3_vector *signature, const mf3 *pk) {

	f3_vector x = f3_vector_new(N - K);
	f3_vector a = f3_vector_new(K);
	f3_vector hash = f3_vector_new(N - K);
	hash_message(&hash, message, mlen, salt, SALT_SIZE);

	/*clock_t start, end;
	 double cpu_time_used;*/

	size_t weight = f3_vector_weight(signature);

	int cmp = 0;
	if (weight != W)
		goto bad_signature;

	memcpy(a.r0, signature->r0 + 45, 88 * sizeof(wave_word));
	memcpy(a.r1, signature->r1 + 45, 88 * sizeof(wave_word));
	shift7(a.r0, a.alloc);
	shift7(a.r1, a.alloc);

	mf3_mv_mul_v(pk, &a, &x);

	f3_vector_sum_inplace_avx(&x, signature);

	size_t nr_words = (N - K) / WORD_LENGTH;

	cmp = memcmp(x.r0, hash.r0, sizeof(wave_word) * nr_words);
	cmp |= memcmp(x.r1, hash.r1, sizeof(wave_word) * nr_words);

	cmp ^= 1;
	bad_signature: f3_vector_free(&x);
	f3_vector_free(&a);
	f3_vector_free(&hash);
	return cmp;

}

/**
 * hash:
 * truncated_signature:
 * pk: we assume that PK is already H^t
 */

int STverify(const uint8_t *salt, const uint8_t *message, const size_t mlen,
		const f3_vector *truncated_signature, unsigned char *sig_compressed,
		size_t size_out, const f3_vector *u, mf3 *pk) {
	int cmp = 0;

	//f3_vector signature_out = f3_vector_new(K);
	//decode_sig(sig_compressed, size_out, &signature_out);

	f3_vector hash = f3_vector_new(N - K);
	hash_message(&hash, message, mlen, salt, SALT_SIZE);

	f3_vector v = f3_vector_new(N - K);
	f3_vector_sub_avx(u, &hash, &v);

	for (int i = 0; i < K; i++) {

		if (f3_vector_coeff(truncated_signature, i) == 2) {
			f3_vector_sum_inplace_avx(&v, &pk->row[i]);
		} else if (f3_vector_coeff(truncated_signature, i) == 0) {
			f3_vector_sub_inplace_avx(&v, &pk->row[i]);

		}
	}

	size_t weight_t = f3_vector_weight(truncated_signature);
	v.r0[v.alloc - 1] = v.r0[v.alloc - 1] & 0x7F;
	v.r1[v.alloc - 1] = v.r1[v.alloc - 1] & 0x7F;
	size_t weight_v = f3_vector_weight(&v);
	/*printf("ST:weight_t: %zu\n", weight_t);
	 printf("ST:weight_v: %zu\n", weight_v);*/

	if ((weight_t + weight_v) == W)
		cmp = 1;

	f3_vector_free(&v);
	f3_vector_free(&hash);

	return cmp;

}

int Nverify(const uint8_t *salt, const uint8_t *message, const size_t mlen,
		const f3_vector *truncated_signature, mf3 *pk) {

	int cmp = 0;
	f3_vector hash = f3_vector_new(N - K);
	hash_message(&hash, message, mlen, salt, SALT_SIZE);

	for (int i = 0; i < K / 2; i++) {
		int c_i, c_ip;
		c_i = (f3_vector_coeff(truncated_signature, 2 * i)
				+ f3_vector_coeff(truncated_signature, (2 * i) + 1)) % 3;
		c_ip =
				((((int) f3_vector_coeff(truncated_signature, 2 * i))
						- (((int) f3_vector_coeff(truncated_signature,
								2 * i + 1)))) + 3) % 3;

		if (c_i == 1) {
			f3_vector_sum_inplace_avx(&hash, &pk->row[2 * i]);
		} else if (c_i == 2) {
			f3_vector_sub_inplace_avx(&hash, &pk->row[2 * i]);
		}
		if (c_ip == 1) {
			f3_vector_sum_inplace_avx(&hash, &pk->row[2 * i + 1]);
		} else if (c_ip == 2) {
			f3_vector_sub_inplace_avx(&hash, &pk->row[2 * i + 1]);
		}

	}

	if (f3_vector_coeff(truncated_signature, K - 1) == 1) {
		f3_vector_sum_inplace_avx(&hash, &pk->row[K - 1]);
	} else if (f3_vector_coeff(truncated_signature, K - 1) == 2) {
		f3_vector_sub_inplace_avx(&hash, &pk->row[K - 1]);
	}
	hash.r0[hash.alloc - 1] = hash.r0[hash.alloc - 1] & 0x7F;
	hash.r1[hash.alloc - 1] = hash.r1[hash.alloc - 1] & 0x7F;

	size_t weight_t = f3_vector_weight(truncated_signature);
	size_t weight_v = f3_vector_weight(&hash);

	if ((weight_t + weight_v) == W)
		cmp = 1;

	f3_vector_free(&hash);

	return cmp;
}

