//	daxpy1(n, a, x, y)
//	do i=1,n
//	  y(i) = y(i) + a * x(i)
//	end do
//
//	modification of daxpy, using pipelined loads of x values
//	reads two elements beyond end of array, assumes strides of 1
//	based on code by Dave Scott, but adapted for use with
//	matrix diagonaliser.
//      bill purvis dl.
//
//
//	r16	N
//	r18	x pointer
//	r19	y pointer
//	r20	loop counter
//	r30	-1 constant for bla
//
//	f2	value of A
//	f4	x1, x5
//	f6	x2, x6
//	f8	x3, x7
//	f10	x4, x8
//	f12	products
//	f14
//	f16	y1
//	f18	y2
//	f20	y3
//	f22	y4
//	f24	y5
//	f26	y6
//	f28	y7
//	f30	y8
//
	.globl	_daxpy1_
	.align	8
_daxpy1_:
	ld.l	r0(r16),r16		// get n
	adds	-1,r0,r30		// -1 for bla
	adds	-8,r19,r19		// correct y pointer
	shr	3,r16,r20		// loop count = n/8
	bte	r20,r0,finale		// bypass main loop for N < 8
	nop				// align for dual mode
	fst.q	f0,-64(sp)++		// save regs
	fst.q	f4,16(sp)
	fst.q	f8,32(sp)
	fst.q	f12,48(sp)
      d.pfiadd.dd f0,f0,f0		// initiate dual mode
	pfld.d	r0(r17),f0		// discard, load A
      d.fnop
	adds	-1,r20,r20		// adjust loop count for bla
      d.fnop
	pfld.d	r0(r18),f0		// discard, 	read x1
      d.fnop
	pfld.d	8(r18)++,f0		// discard,	read x2
      d.fnop
	bla	r30,r20,loop		// prime LCC
      d.fnop
	pfld.d	8(r18)++,f2		// load A,	read x3
loop: d.fnop
	fld.d	8(r19),f16		// read y1
      d.fnop
	pfld.d	8(r18)++,f4		// load x1,	read x4
      d.fnop
	fld.d	16(r19),f18		// read y2
      d.fnop
	pfld.d	8(r18)++,f6		// load x2,	read x5
      d.pfmul.dd f2,f4,f0		// -		start a*x1
	fld.d	24(r19),f20		// read y3
      d.fnop
	pfld.d	8(r18)++,f8		// load x3,	read x6
      d.pfmul.dd f2,f6,f0		// -		start a*x2
	fld.d	32(r19),f22		// read y4
      d.fnop
	pfld.d	8(r18)++,f10		// load x4,	read x7
      d.pfmul.dd f2,f8,f12		// store a*x1,	start a*x3
	fld.d	40(r19),f24		// read y5
      d.pfadd.dd f16,f12,f0		// -		start (a*x1)+y1
	pfld.d	8(r18)++,f4		// load x5,	read x8
      d.pfmul.dd f2,f10,f12		// store a*x2,	start a*x4
	fld.d	48(r19),f26		// read y6
      d.pfadd.dd f18,f12,f0		// -		start (a*x2)+y2
	pfld.d	8(r18)++,f6		// load x6,	read x1'
      d.pfmul.dd f2,f4,f12		// store a*x3,	start a*x5
	fld.d	56(r19),f28      	// read y7
      d.pfadd.dd f20,f12,f0		// -		start (a*x3)+y3
	pfld.d	8(r18)++,f8		// load x7,	read x2'
      d.pfmul.dd f2,f6,f12		// store a*x4,  start a*x6
	fld.d	64(r19),f30		// read y8
      d.pfadd.dd f22,f12,f16		// store y1	start (a*x4)+y4
	pfld.d	8(r18)++,f10		// load x8,	read x3'
      d.pfmul.dd f2,f8,f12		// store a*x5,	start a*x7
	fst.d	f16,8(r19)++		// save y1
      d.pfadd.dd f24,f12,f18		// store y2	start (a*x5)+y5
	nop
      d.pfmul.dd f2,f10,f12		// store a*x6	start a*x8
	nop
      d.pfadd.dd f26,f12,f20		// store y3	start (a*x6)+y6
	nop
      d.pfmul.dd f0,f0,f12		// store a*x7	flush
	fst.d	f18,8(r19)++		// save y2
      d.pfadd.dd f28,f12,f22		// store y4,	start (a*x7)+y7
	fst.d	f20,8(r19)++		// save y3
      d.pfmul.dd f0,f0,f12		// store a*x8	flush
	fst.d	f22,8(r19)++		// save y4
      d.pfadd.dd f30,f12,f24		// store y5,	start (a*x8)+y8
	fst.d	f24,8(r19)++		// save y5
      d.pfadd.dd f0,f0,f26		// store y6,	flush
	fst.d	f26,8(r19)++		// save y6
      d.pfadd.dd f0,f0,f28		// store y7,	flush
        fst.d	f28,8(r19)++		// save y7
      d.pfadd.dd f0,f0,f30		// store y8,	flush
	bla	r30,r20,loop		// loop control
      d.fnop
	fst.d	f30,8(r19)++		// save y8
//
//	exit dual mode, restore fp regs
//
	fnop
	fld.q	0(sp),f0
	fnop
	fld.q	16(sp),f4
	fld.q	32(sp),f8
	fld.q	48(sp),f12
	adds	64,sp,sp
	adds	-16,r18,r18		// correct overrun for remainder
finale:
	and	7,r16,r16		// mask off to get remainder
	bc	exit			// skip if nothing left	
	adds	r30,r16,r16		// pre-adjust counter
	adds	-8,r18,r18		// adjust x pointer
	bla	r30,r16,loop2
	fld.d	r0(r17),f16		// re-load a
loop2:
	fld.d	8(r18)++,f18		// load x value
	fld.d	8(r19),f20		// load y value
	fmul.dd	f16,f18,f18		// a * x(i)
	fadd.dd	f18,f20,f20		// + y(i)
	bla	r30,r16,loop2
	fst.d	f20,8(r19)++		// ...storing updated value
//
exit:
	bri	r1
	nop
