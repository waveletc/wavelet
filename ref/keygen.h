/*
 * keygen.h
 *
 *  Created on: Jan 29, 2021
 *      Author: vader
 */

#ifndef KEYGEN_H_
#define KEYGEN_H_

#include <stdio.h>
#include <stdint.h>

#include "f3.h"
#include "vf3.h"
#include "mf3.h"

#include "prng.h"
#include "randombytes.h"

#include "params.h"
#include "definitions.h"
#include "enum.h"


void phiUV(int8_t *e, int8_t *eU, int8_t *eV, wave_sk_t sk);

int keygen(wave_sk_t *sk, wave_pk_t *pk);

void wave_sk_clear(wave_sk_t *sk);

void wave_pk_clear(wave_pk_t pk);

#endif /* KEYGEN_H_ */
