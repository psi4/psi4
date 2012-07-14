/*-
 * Copyright (c) 2012 Ilya Kaliman
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

static const char *potential_files[] = {
	ABS_TOP_SRCDIR "/fraglib/h2o.efp",
	ABS_TOP_SRCDIR "/fraglib/nh3.efp",
	ABS_TOP_SRCDIR "/fraglib/ch3oh.efp",
	NULL
};

static const char *fragname[] = {
	"H2O_L",
	"NH3_L",
	"NH3_L",
	"NH3_L",
	"CH3OH_L",
	"H2O_L",
	"H2O_L",
	"CH3OH_L",
	"H2O_L",
	 NULL
};

static const double geometry[] = {
	BOHR(-3.394), BOHR(-1.900), BOHR(-3.700), /* H2O */
	BOHR(-3.524), BOHR(-1.089), BOHR(-3.147),
	BOHR(-2.544), BOHR(-2.340), BOHR(-3.445),
	BOHR(-5.515), BOHR( 1.083), BOHR( 0.968), /* NH3 */
	BOHR(-5.161), BOHR( 0.130), BOHR( 0.813),
	BOHR(-4.833), BOHR( 1.766), BOHR( 0.609),
	BOHR( 1.848), BOHR( 0.114), BOHR( 0.130), /* NH3 */
	BOHR( 1.966), BOHR( 0.674), BOHR(-0.726),
	BOHR( 0.909), BOHR( 0.273), BOHR( 0.517),
	BOHR(-1.111), BOHR(-0.084), BOHR(-4.017), /* NH3 */
	BOHR(-1.941), BOHR( 0.488), BOHR(-3.813),
	BOHR(-0.292), BOHR( 0.525), BOHR(-4.138),
	BOHR(-2.056), BOHR( 0.767), BOHR(-0.301), /* CH3OH */
	BOHR(-2.999), BOHR(-0.274), BOHR(-0.551),
	BOHR(-1.201), BOHR( 0.360), BOHR( 0.258),
	BOHR(-0.126), BOHR(-2.228), BOHR(-0.815), /* H2O */
	BOHR( 0.310), BOHR(-2.476), BOHR( 0.037),
	BOHR( 0.053), BOHR(-1.277), BOHR(-1.011),
	BOHR(-1.850), BOHR( 1.697), BOHR( 3.172), /* H2O */
	BOHR(-1.050), BOHR( 1.592), BOHR( 2.599),
	BOHR(-2.666), BOHR( 1.643), BOHR( 2.614),
	BOHR( 1.275), BOHR(-2.447), BOHR(-4.673), /* CH3OH */
	BOHR( 0.709), BOHR(-3.191), BOHR(-3.592),
	BOHR( 2.213), BOHR(-1.978), BOHR(-4.343),
	BOHR(-5.773), BOHR(-1.738), BOHR(-0.926), /* H2O */
	BOHR(-5.017), BOHR(-1.960), BOHR(-1.522),
	BOHR(-5.469), BOHR(-1.766), BOHR( 0.014)
};
