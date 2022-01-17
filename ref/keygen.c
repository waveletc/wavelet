#include <stdlib.h>
#include <stdio.h>
#include "keygen.h"

void gen_coeffUV(prng_t *PRNG, coeff_t *coeff) {
	int i;

	for (i = 0; i < N2; i++) {
		coeff->a[i] = ((1 - (2 * rnd_bit(PRNG)) + 3) % 3);
		coeff->c[i] = ((1 - (2 * rnd_bit(PRNG)) + 3) % 3);

		coeff->b[i] = rnd_trit(PRNG);
		coeff->d[i] = f3_mul(coeff->a[i],
				f3_add(1, f3_mul(coeff->b[i], coeff->c[i])));
		coeff->x[i] = f3_add(f3_mul(coeff->a[i], coeff->b[i]),
				f3_mul(coeff->c[i], coeff->d[i]));
	}

}

mf3* full_parity_check_matrix(mf3 *HU, mf3 *HV, coeff_t *coeff) {
	int i, j;
	mf3 *Hsec = mf3_new(N - K, N);

	for (i = 0; i < N2 - KU; i++) {
		for (j = 0; j < N2; j++) {
			mf3_setcoeff(Hsec, i, j, f3_mul(mf3_coeff(HU, i, j), coeff->d[j]));
			mf3_setcoeff(Hsec, i, j + N2,
					f3_neg(f3_mul(mf3_coeff(HU, i, j), coeff->b[j])));

		}
	}
	for (i = 0; i < N2 - KV; i++) {
		for (j = 0; j < N2; j++) {

			mf3_setcoeff(Hsec, i + N2 - KU, j,
					f3_neg(f3_mul(mf3_coeff(HV, i, j), coeff->c[j])));

			mf3_setcoeff(Hsec, i + N2 - KU, j + N2,
					f3_mul(mf3_coeff(HV, i, j), coeff->a[j]));

		}
	}

	return Hsec;
}

int keygen(wave_sk_t *sk, wave_pk_t *pk) {

	uint8_t key_seed = 2;
	randombytes(&key_seed, 1);
	prng_t *PRNG = prng_init(key_seed);
	randperm(sk->perm, N);

	sk->HU = mf3_rand(N2 - KU, N2, PRNG);

	sk->HV = mf3_rand(N2 - KV, N2, PRNG);

	gen_coeffUV(PRNG, &sk->coeff);

	mf3 *H = full_parity_check_matrix(sk->HU, sk->HV, &sk->coeff);

	mf3 *Hsec = mf3_copy(H);

	mf3_gauss_elim(H, sk->perm);

	mf3 *R = mf3_new(N - K, K);

	for (int i = 0; i < N - K; i++) {
		for (int j = 0; j < K; j++) {

			f3_vector_setcoeff(&R->row[i], j,
					mf3_coeff(H, i, sk->perm[N - K + j]));
		}
	}

	//TODO: recompute it at sign?
	sk->Sinv = mf3_new(N - K, N - K);
	for (int i = 0; i < (N - K); i++) {
		for (int j = 0; j < (N - K); j++) {
			f3_vector_setcoeff(&sk->Sinv->row[i], j,
					mf3_coeff(Hsec, i, sk->perm[j]));

		}
	}

	*pk = R;

	mf3_free(H);
	mf3_free(Hsec);
	prng_clear(PRNG);

	return SUCCESS;
}

void wave_sk_clear(wave_sk_t *sk) {
	mf3_free(sk->HU);
	mf3_free(sk->HV);
	mf3_free(sk->Sinv);
}

void wave_pk_clear(wave_pk_t pk) {
	mf3_free(pk);
}

