	.file	"/home/d3g681/tcgmsg/ipcv4.0/"
	.file	"copyall.c"
	.vstamp 7
# KSR1 ccom  -OLM -X28 -X92 -X115 -X151 -X153 -X155 -X156 -X157 -X158 -X159
# 	 -X172 -X187
# ccom: version 1.1.1. built Sun Dec 26 22:03:57 1993.

	.text

	.data
	.def	copyto$TXT;	.val	copyto$TXT;	.scl	2;	.endef

	.text
	.def	copyto;	.val	copyto;	.scl	2;	.type	513;	.endef
copyto$TXT:
   finop				;  cxnop
   finop				;  cxnop
	.def	.bf;	.val	.;	.scl	101;	.line	11;	.endef
   mov8_8	%i3, %i9		;  ssub8.ntr	0, %sp, 128, %sp
   itstle8	128, %i4		;  movb8_8	%i2, %c8
   add8.ntr	10, %i31, %i31		;  st8		%i13, 80(%sp)
   mov8_8	%i4, %i13		;  st8		%cp, 112(%sp)
   finop				;  st8		%fp, 120(%sp)
   finop				;  mov8_8	%c10, %cp
   finop				;  sadd8.ntr	0, %sp, 128, %fp
   finop				;  bcc.qn	@citst, .L2
   finop				;  st8		%c14, 104(%sp)
   finop				;  st8		%i12, 88(%sp)
   clrh8	7, %i9, %i5		;  movb8_8	%c8, %i1
   sub8.ntr	128, %i5, %i2		;  cxnop
   add8.ntr	7, %i31, %i31		;  cxnop
   sub8.ntr	%i9, %i1, %i1		;  cxnop
   clrh8	3, %i1, %i1		;  cxnop
   itsteq8	0, %i1			;  cxnop
   itstge8	0, %i5			;  bcs.qt	@citst, .L10
.L2:
   mov8_8	%i13, %i4		;  ld8		16(%cp), %c6
   mov8_8	%i9, %i2		;  movb8_8	%c8, %i3
   add8.ntr	4, %i31, %i31		;  ld8		8(%cp), %c10
   finop				;  jsr		%c14, 16(%c6)
   finop				;  cxnop
   finop				;  cxnop
   movi8	3, %i0			;  movi8	0, %c8
	.ln	7, .-32	# 17

   add8.ntr	8, %i31, %i31		;  ld8		104(%sp), %c14
   finop				;  ld8		112(%sp), %cp
   finop				;  ld8		120(%sp), %fp
   finop				;  ld8		88(%sp), %i12
   finop				;  ld8		80(%sp), %i13
   finop				;  jmp		32(%c14)
   finop				;  sadd8.ntr	0, %sp, 128, %sp
   finop				;  cxnop
.L10:
   selsc8	%i5, %i2, %i5		;  cxnop
   itstle8	%i5, %i13		;  cxnop
   selsc8	%i5, %i13, %i5		;  cxnop
   sub8.ntr	%i13, %i5, %i13		;  cxnop
   itsteq8	0, %i5			;  cxnop
   sub8.ntr	%i5, 1, %i5		;  bcs.qt	@citst, .L8
.L9:
   add8.ntr	1, %i9, %i9		;  movb8_8	%i9, %c4
   itsteq8	0, %i5			;  ld1		0(%c8), %i10
   sub8.ntr	%i5, 1, %i5		;  bcc.qn	@citst, .L9
   add8.ntr	3, %i31, %i31		;  sadd8.ntr	0, %c8, 1, %c8
   finop				;  st1		%i10, -1(%c4)
.L8:
   itstne8	0, %i13			;  movb8_8	%i9, %c7
   ash8.ntr	-7, %i13, %i10		;  bcc.qt	@citst, .L1
   ash8.ntr	7, %i10, %i0		;  movb8_8	%i0, %c5
   add8.ntr	%i9, %i0, %i9		;  mov8_8	%c8, %c6
   finop				;  pcsp.ex.bl  128(%c7)
   finop				;  pcsp.ex.bl  256(%c7)
   finop				;  pcsp.ex.bl  384(%c7)
   itsteq8	0, %i10			;  movb8_8	%i9, %c9
   sub8.ntr	%i13, %i0, %i13		;  sadd8.ntr	0, %c5, %c8, %c8
   sub8.ntr	%i10, 1, %i10		;  bcs.qt	@citst, .L5
.L6:
   finop				;  pcsp.ex.bl   512(%c7)
   itsteq8	0, %i10			;  ld8.ro	0(%c6), %i11
   sub8.ntr	%i10, 1, %i10		;  ld8.ro	8(%c6), %i0
   add8.ntr	33, %i31, %i31		;  ld8.ro	16(%c6), %i1
   finop				;  ld8.ro	24(%c6), %i2
   finop				;  ld8.ro	32(%c6), %i3
   finop				;  ld8.ro	40(%c6), %i4
   finop				;  ld8.ro	48(%c6), %i5
   finop				;  ld8.ro	56(%c6), %i12
   finop				;  sadd8.ntr	0, %c7, 128, %c7
   finop				;  st8		%i11, -128(%c7)
   finop				;  st8		%i0, -120(%c7)
   finop				;  st8		%i1, -112(%c7)
   finop				;  st8		%i2, -104(%c7)
   finop				;  st8		%i3, -96(%c7)
   finop				;  st8		%i4, -88(%c7)
   finop				;  st8		%i5, -80(%c7)
   finop				;  st8		%i12, -72(%c7)
   finop				;  ld8.ro	120(%c6), %i12
   finop				;  ld8.ro	112(%c6), %i5
   finop				;  ld8.ro	104(%c6), %i4
   finop				;  ld8.ro	96(%c6), %i3
   finop				;  ld8.ro	88(%c6), %i2
   finop				;  ld8.ro	80(%c6), %i1
   finop				;  ld8.ro	72(%c6), %i0
   finop				;  ld8.ro	64(%c6), %i11
   finop				;  st8		%i1, -48(%c7)
   finop				;  st8		%i2, -40(%c7)
   finop				;  st8		%i0, -56(%c7)
   finop				;  st8		%i11, -64(%c7)
   finop				;  st8		%i3, -32(%c7)
   finop				;  st8		%i4, -24(%c7)
   finop				;  st8		%i5, -16(%c7)
   finop				;  bcc.qn	@citst, .L6
   finop				;  st8		%i12, -8(%c7)
   finop				;  sadd8.ntr	0, %c6, 128, %c6
#   finop				;  pstsp        0(%c7) 
.L5:
   itsteq8	0, %i13			;  cxnop
   sub8.ntr	%i13, 1, %i4		;  bcs.qt	@citst, .L1
.L4:
   itsteq8	0, %i4			;  ld1		0(%c8), %i11
   sub8.ntr	%i4, 1, %i4		;  sadd8.ntr	0, %c9, 1, %c9
   add8.ntr	3, %i31, %i31		;  bcc.qn	@citst, .L4
   finop				;  sadd8.ntr	0, %c8, 1, %c8
   finop				;  st1		%i11, -1(%c9)
.L1:
   add8.ntr	8, %i31, %i31		;  ld8		104(%sp), %c14
   finop				;  ld8		112(%sp), %cp
   finop				;  ld8		120(%sp), %fp
   finop				;  ld8		88(%sp), %i12
   finop				;  ld8		80(%sp), %i13
   finop				;  jmp		32(%c14)
   finop				;  sadd8.ntr	0, %sp, 128, %sp
   finop				;  cxnop
	.def	.ef;	.val	.;	.scl	101;	.line	93;	.endef
	.def	copyto;	.scl	-1;	.endef

	.data
# nbytes	%i5	local
# npage	%i10	local
# from	%c6	local
# to	%c7	local
# a	%i11	local
# b	%i0	local
# c	%i1	local
# d	%i2	local
# e	%i3	local
# f	%i4	local
# g	%i5	local
# h	%i12	local
# nbytes	%i4	local
# from	%c8	local
# to	%c9	local
	.half	  0x0, 0x0, 0x60003000, 0x5800
.L21:
copyto:	.word	  copyto$TXT
	.word	  memcpy
	.word	  memcpy$TXT

# src	%c8	local
# dest	%i9	local
# n	%i13	local

	.text

	.data
	.def	copyfrom$TXT;	.val	copyfrom$TXT;	.scl	2;	.endef

	.text
	.def	copyfrom;	.val	copyfrom;	.scl	2;	.type	513;	.endef
copyfrom$TXT:
   finop				;  cxnop
   finop				;  cxnop
	.def	.bf;	.val	.;	.scl	101;	.line	112;	.endef
   itstle8	128, %i4		;  ssub8.ntr	0, %sp, 128, %sp
   add8.ntr	10, %i31, %i31		;  movb8_8	%i3, %c8
   finop				;  st8		%i13, 80(%sp)
   mov8_8	%i4, %i13		;  st8		%cp, 112(%sp)
   finop				;  st8		%fp, 120(%sp)
   finop				;  mov8_8	%c10, %cp
   finop				;  sadd8.ntr	0, %sp, 128, %fp
   finop				;  bcc.qn	@citst, .L25
   finop				;  st8		%c14, 104(%sp)
   finop				;  st8		%i12, 88(%sp)
   clrh8	7, %i2, %i5		;  movb8_8	%c8, %i0
   sub8.ntr	128, %i5, %i1		;  cxnop
   add8.ntr	7, %i31, %i31		;  cxnop
   sub8.ntr	%i0, %i2, %i0		;  cxnop
   clrh8	3, %i0, %i0		;  cxnop
   itsteq8	0, %i0			;  cxnop
   itstge8	0, %i5			;  bcs.qt	@citst, .L33
.L25:
   mov8_8	%i2, %i3		;  ld8		16(%cp), %c6
   mov8_8	%i13, %i4		;  movb8_8	%c8, %i2
   add8.ntr	4, %i31, %i31		;  ld8		8(%cp), %c10
   finop				;  jsr		%c14, 16(%c6)
   finop				;  cxnop
   finop				;  cxnop
   movi8	3, %i0			;  movi8	0, %c8
	.ln	7, .-32	# 118

   add8.ntr	8, %i31, %i31		;  ld8		104(%sp), %c14
   finop				;  ld8		112(%sp), %cp
   finop				;  ld8		120(%sp), %fp
   finop				;  ld8		88(%sp), %i12
   finop				;  ld8		80(%sp), %i13
   finop				;  jmp		32(%c14)
   finop				;  sadd8.ntr	0, %sp, 128, %sp
   finop				;  cxnop
.L33:
   selsc8	%i5, %i1, %i5		;  cxnop
   itstle8	%i5, %i13		;  cxnop
   selsc8	%i5, %i13, %i5		;  cxnop
   sub8.ntr	%i13, %i5, %i13		;  cxnop
   itsteq8	0, %i5			;  cxnop
   sub8.ntr	%i5, 1, %i5		;  bcs.qt	@citst, .L31
.L32:
   add8.ntr	1, %i2, %i2		;  movb8_8	%i2, %c4
   itsteq8	0, %i5			;  sadd8.ntr	0, %c8, 1, %c8
   sub8.ntr	%i5, 1, %i5		;  cxnop
   add8.ntr	5, %i31, %i31		;  ld1		-1(%c4), %i9
   finop				;  bcc.qn	@citst, .L32
   finop				;  cxnop
   finop				;  st1		%i9, -1(%c8)
.L31:
   itstne8	0, %i13			;  movb8_8	%i2, %c6
   finop				;  pcsp.ro.bl  128(%c6)
   finop				;  pcsp.ro.bl  256(%c6)
   finop				;  pcsp.ro.bl  384(%c6)
   ash8.ntr	-7, %i13, %i9		;  bcc.qt	@citst, .L24
   ash8.ntr	7, %i9, %i11		;  movb8_8	%i11, %c5
   add8.ntr	%i2, %i11, %i2		;  mov8_8	%c8, %c7
   itsteq8	0, %i9			;  movb8_8	%i2, %c9
   sub8.ntr	%i13, %i11, %i13	;  sadd8.ntr	0, %c5, %c8, %c8
   sub8.ntr	%i9, 1, %i9		;  bcs.qt	@citst, .L28
.L29:
   finop				;  pcsp.ro.bl  512(%c6)
   itsteq8	0, %i9			;  ld8.ro	0(%c6), %i10
   sub8.ntr	%i9, 1, %i9		;  ld8.ro	8(%c6), %i11
   add8.ntr	33, %i31, %i31		;  ld8.ro	16(%c6), %i0
   finop				;  ld8.ro	24(%c6), %i1
   finop				;  ld8.ro	32(%c6), %i3
   finop				;  ld8.ro	40(%c6), %i4
   finop				;  ld8.ro	48(%c6), %i5
   finop				;  ld8.ro	56(%c6), %i12
   finop				;  sadd8.ntr	0, %c7, 128, %c7
   finop				;  st8		%i10, -128(%c7)
   finop				;  st8		%i11, -120(%c7)
   finop				;  st8		%i0, -112(%c7)
   finop				;  st8		%i1, -104(%c7)
   finop				;  st8		%i3, -96(%c7)
   finop				;  st8		%i4, -88(%c7)
   finop				;  st8		%i5, -80(%c7)
   finop				;  st8		%i12, -72(%c7)
   finop				;  ld8.ro	120(%c6), %i12
   finop				;  ld8.ro	112(%c6), %i5
   finop				;  ld8.ro	104(%c6), %i4
   finop				;  ld8.ro	96(%c6), %i3
   finop				;  ld8.ro	88(%c6), %i1
   finop				;  ld8.ro	80(%c6), %i0
   finop				;  ld8.ro	72(%c6), %i11
   finop				;  ld8.ro	64(%c6), %i10
   finop				;  st8		%i0, -48(%c7)
   finop				;  st8		%i1, -40(%c7)
   finop				;  st8		%i11, -56(%c7)
   finop				;  st8		%i10, -64(%c7)
   finop				;  st8		%i3, -32(%c7)
   finop				;  st8		%i4, -24(%c7)
   finop				;  st8		%i5, -16(%c7)
   finop				;  bcc.qn	@citst, .L29
   finop				;  sadd8.ntr	0, %c6, 128, %c6
   finop				;  st8		%i12, -8(%c7)
.L28:
   itsteq8	0, %i13			;  cxnop
   sub8.ntr	%i13, 1, %i4		;  bcs.qt	@citst, .L24
.L27:
   itsteq8	0, %i4			;  ld1		0(%c9), %i10
   sub8.ntr	%i4, 1, %i4		;  sadd8.ntr	0, %c8, 1, %c8
   add8.ntr	3, %i31, %i31		;  bcc.qn	@citst, .L27
   finop				;  sadd8.ntr	0, %c9, 1, %c9
   finop				;  st1		%i10, -1(%c8)
.L24:
   add8.ntr	8, %i31, %i31		;  ld8		104(%sp), %c14
   finop				;  ld8		112(%sp), %cp
   finop				;  ld8		120(%sp), %fp
   finop				;  ld8		88(%sp), %i12
   finop				;  ld8		80(%sp), %i13
   finop				;  jmp		32(%c14)
   finop				;  sadd8.ntr	0, %sp, 128, %sp
   finop				;  cxnop
	.def	.ef;	.val	.;	.scl	101;	.line	93;	.endef
	.def	copyfrom;	.scl	-1;	.endef

	.data
# nbytes	%i5	local
# npage	%i9	local
# from	%c6	local
# to	%c7	local
# a	%i10	local
# b	%i11	local
# c	%i0	local
# d	%i1	local
# e	%i3	local
# f	%i4	local
# g	%i5	local
# h	%i12	local
# nbytes	%i4	local
# from	%c9	local
# to	%c8	local
	.half	  0x0, 0x0, 0x60003000, 0x5800
.L44:
copyfrom:	.word	  copyfrom$TXT
	.word	  memcpy
	.word	  memcpy$TXT

# src	%i2	local
# dest	%c8	local
# n	%i13	local

	.text

	.data

	.align  	128
.L47:
	.globl  	copyfrom
	.globl  	copyfrom$TXT
	.globl  	copyto
	.globl  	copyto$TXT

	.text
