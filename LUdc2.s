	# LU decomposition of a matrix by D.L.Buunk
	# for x86-64 processors with SSE2

	.text
	.global LUdecomp2
LUdecomp2:
	# save registers
	pushq	%rbx
	pushq	%rcx
	pushq	%rdx

	pushq	%r8
	pushq	%r9
	pushq	%r10
	pushq	%r11

	pushq	%rsi
	pushq	%rdi
	pushq	%rbp

	# set up variables in their proper registers
	movq	%rdx,%rbp	# * permutations
	movq	%rcx,%rdx	# dimensions

	# set up permutation array (zero-based)
	xorq	%rax,%rax
	movq	%rdx,%rcx	# dimensions
	movq	%rbp,%rbx	# * permutations
setup_perm:
	movq	%rax,(%rbx)
	addq	$0x4,%rbx
	incq	%rax
	loop	setup_perm

	# set up pointer array a
	movq	%rdx,%rcx	# dimensions
	movq	%rdi,%rbx	# * matrix a
	movq	%rsi,%r8	# ** matrix a
	movq	%rdx,%rax	# dimensions
	# take care of alignment
	decq	%rax
	shrq	$1,%rax
	incq	%rax
	shlq	$3,%rax
setup_a:
	movq	%rbx,(%rsi)
	addq	$0x8,%rsi
	addq	%rax,%rbx
	loop	setup_a

	# set up pointer array b
	movq	%rdx,%rcx	# dimensions
	movq	%rsi,%r9	# ** matrix bT
	movq	%rdx,%rax	# dimensions
	shlq	$2,%rax
	addq	%rsi,%rax	# * matrix bT
	xorq	%rbx,%rbx	# size of a row of bT
setup_b:
	addq	$0x10,%rbx

	movq	%rax,(%rsi)
	addq	$0x8,%rsi
	addq	%rbx,%rax
	decq	%rcx
	jrcxz	setup_end

	movq	%rax,(%rsi)
	addq	$0x8,%rsi
	addq	%rbx,%rax
	decq	%rcx
	jrcxz	setup_end

	movq	%rax,(%rsi)
	addq	$0x8,%rsi
	addq	%rbx,%rax
	decq	%rcx
	jrcxz	setup_end

	movq	%rax,(%rsi)
	addq	$0x8,%rsi
	addq	%rbx,%rax
	loop	setup_b
setup_end:
	# okay, setup is done, now we can setup the outer loop
	xorq	%rax,%rax
	xorq	%rbx,%rbx
	# than, jump to pivoting, as we can skip the beta loop for the first col
	jmp	pivoting
outer_loop:

	# inside it there is the beta loop, the pivoting, and the alpha loop
	xorq	%rbx,%rbx
beta_loop:

	# get the pointers for the dot-product
	movq	%rbx,%rcx
	shlq	$3,%rcx
	addq	%r8,%rcx
	movq	(%rcx),%rdi

	movq	%rax,%rcx
	shlq	$3,%rcx
	addq	%r9,%rcx
	movq	(%rcx),%rsi

	# clear xmm0
	xorps	%xmm0,%xmm0

	# inside the beta loop, there is an inner loop for the dot-product
	movq	%rbx,%rcx
inner_loop_a:
	cmpq	$4,%rcx
	jb	inner_a_almost

	# multiply and add 4 floats at the same time
	movaps	(%rdi),%xmm1
	movaps	(%rsi),%xmm2
	mulps	%xmm2,%xmm1
	addq	$0x10,%rdi
	addq	$0x10,%rsi
	addps	%xmm1,%xmm0

	subq	$4,%rcx
	jmp	inner_loop_a

inner_a_almost:
	jrcxz	inner_a_0

	# 3/2/1 left
	xorps	%xmm1,%xmm1
	xorps	%xmm2,%xmm2
	movss	(%rdi),%xmm1
	movss	(%rsi),%xmm2
	decq	%rcx
	jrcxz	inner_a_end

	# 2/1 left
	pslldq	$4,%xmm1
	pslldq	$4,%xmm2
	addq	$0x4,%rdi
	addq	$0x4,%rsi
	movss	(%rdi),%xmm1
	movss	(%rsi),%xmm2
	decq	%rcx
	jrcxz	inner_a_end

	# 1 left
	pslldq	$4,%xmm1
	pslldq	$4,%xmm2
	addq	$0x04,%rdi
	addq	$0x04,%rsi
	movss	(%rdi),%xmm1
	movss	(%rsi),%xmm2

inner_a_end:
	addq	$0x04,%rsi
	mulps	%xmm2,%xmm1
	addps	%xmm1,%xmm0

inner_a_0:
	# now we only have to add up the 4 values to 1
	movaps	%xmm0,%xmm1
	psrldq	$8,%xmm1
	addps	%xmm1,%xmm0
	movaps	%xmm0,%xmm1
	psrldq	$4,%xmm1
	addss	%xmm1,%xmm0

	# now, get the "a" value and substract.
	# calculate the pointer first, we can use %rcx
	movq	%rbx,%rcx
	shlq	$3,%rcx
	addq	%r8,%rcx
	movq	(%rcx),%rdi
	movq	%rax,%rcx
	shlq	$2,%rcx
	addq	%rcx,%rdi
	movss	(%rdi),%xmm1
	subss	%xmm0,%xmm1
	# store in "a" matrix
	movss	%xmm1,(%rdi)
	# store in "bT" matrix
	movss	%xmm1,(%rsi)

	# end of beta loop
	incq	%rbx
	cmpq	%rbx,%rax
	jne	beta_loop

pivoting:
	# pivoting runs from k=i to k=N
	# r10 is current, r11 is best, rcx is counter
	# xmm3[0:31] is the best val
	movq	$-1,%r11
	movq	%rbx,%r10
	xorps	%xmm3,%xmm3

pivot_loop:
	# get the pointers to the matrices

	movq	%rax,%rcx
	shlq	$3,%rcx
	addq	%r9,%rcx
	movq	(%rcx),%rsi

	movq	%r10,%rcx
	shlq	$3,%rcx
	addq	%r8,%rcx
	movq	(%rcx),%rdi

	# clear xmm0
	xorps	%xmm0,%xmm0

	# inner loop
	movq	%rdx,%rcx
	subq	%rbx,%rcx
inner_loop_p:
	cmpq	$4,%rcx
	jb	inner_p_almost

	# mul and add 4 floats
	movaps	(%rdi),%xmm1
	movaps	(%rsi),%xmm2
	mulps	%xmm2,%xmm1
	addq	$0x10,%rdi
	addq	$0x10,%rsi
	addps	%xmm1,%xmm0

	subq	$4,%rcx
	jmp	inner_loop_p

inner_p_almost:
	jrcxz	inner_p_0

	# 3/2/1 left
	xorps	%xmm1,%xmm1
	xorps	%xmm2,%xmm2
	movss	(%rdi),%xmm1
	movss	(%rsi),%xmm2
	decq	%rcx
	jrcxz	inner_p_end

	# 2/1 left
	pslldq	$4,%xmm1
	pslldq	$4,%xmm2
	addq	$0x04,%rdi
	addq	$0x04,%rsi
	movss	(%rdi),%xmm1
	movss	(%rsi),%xmm2
	decq	%rcx
	jrcxz	inner_p_end

	# 1 left
	pslldq	$4,%xmm1
	pslldq	$4,%xmm2
	addq	$0x04,%rdi
	addq	$0x04,%rsi
	movss	(%rdi),%xmm1
	movss	(%rsi),%xmm2

inner_p_end:
	mulps	%xmm2,%xmm1
	addps	%xmm1,%xmm0

inner_p_0:
	# final add
	movaps	%xmm0,%xmm1
	psrldq	$8,%xmm1
	addps	%xmm1,%xmm0
	movaps	%xmm0,%xmm1
	psrldq	$4,%xmm1
	addss	%xmm1,%xmm0

	# get the "a" val and substract
	# find the pointer first
	movq	%r10,%rcx
	shlq	$3,%rcx
	addq	%r8,%rcx
	movq	(%rcx),%rdi
	movq	%rax,%rcx
	shlq	$2,%rcx
	addq	%rcx,%rdi
	movss	(%rdi),%xmm1
	subss	%xmm0,%xmm1

	# absolute value
	xorps	%xmm0,%xmm0
	movq	$0x7FFFFFFF,%rcx
	pushq	%rcx
	movss	(%rsp),%xmm0
	popq	%rcx
	pand	%xmm0,%xmm1

	# get maximal value
	movss	%xmm3,%xmm0
	subss	%xmm1,%xmm0
	movmskps %xmm0,%rcx
	andq	$1,%rcx
	cmpq	$1,%rcx
	jne	end_pivot

	# store max val and index
	movss	%xmm1,%xmm3
	movq	%r10,%r11

end_pivot:
	# end of pivot loop
	incq	%r10
	cmpq	%r10,%rdx
	jne	pivot_loop

	# swap pointers
	movq	%rbx,%rdi
	movq	%r11,%rsi
	shlq	$3,%rdi
	shlq	$3,%rsi
	addq	%r8,%rdi
	addq	%r8,%rsi
	movq	(%rsi),%rcx
	xchgq	(%rdi),%rcx
	movq	%rcx,(%rsi)

	# swap indices in perm-array
	movq	%rbx,%rdi
	movq	%r11,%rsi
	shlq	$2,%rdi
	shlq	$2,%rsi
	addq	%rbp,%rdi
	addq	%rbp,%rsi
	movl	(%rsi),%ecx
	xchgl	(%rdi),%ecx
	movl	%ecx,(%rsi)

	# store the value to both arrays
	# as always, first get the pointers
	movq	%rbx,%rcx
	shlq	$3,%rcx
	addq	%r8,%rcx
	movq	(%rcx),%rdi
	movq	%rax,%rcx
	shlq	$2,%rcx
	addq	%rcx,%rdi
	movss	%xmm3,(%rdi)
	movq	%rax,%rcx
	shlq	$3,%rcx
	addq	%r8,%rcx
	movq	(%rcx),%rsi
	movq	%rbx,%rcx
	shlq	$2,%rcx
	addq	%rcx,%rsi
	movss	%xmm3,(%rsi)

alpha_loop:

	# get the pointers for the dot-product
	movq	%rbx,%rcx
	shlq	$3,%rcx
	addq	%r8,%rcx
	movq	(%rcx),%rdi

	movq	%rax,%rcx
	shlq	$3,%rcx
	addq	%r9,%rcx
	movq	(%rcx),%rsi

	# clear xmm0
	xorps	%xmm0,%xmm0

	movq	%rax,%rcx
inner_loop_b:
	cmpq	$4,%rcx
	jb	inner_b_almost

	# multiply and add
	movaps	(%rdi),%xmm1
	movaps	(%rsi),%xmm2
	mulps	%xmm2,%xmm1
	addq	$0x10,%rdi
	addq	$0x10,%rsi
	addps	%xmm1,%xmm0

	subq	$4,%rcx
	jmp	inner_loop_b

inner_b_almost:
	jrcxz	inner_b_0

	# 3/2/1 left
	xorps	%xmm1,%xmm1
	xorps	%xmm2,%xmm2
	movss	(%rdi),%xmm1
	movss	(%rsi),%xmm2
	decq	%rcx
	jrcxz	inner_b_end

	# 2/1 left
	pslldq	$4,%xmm1
	pslldq	$4,%xmm2
	addq	$0x4,%rdi
	addq	$0x4,%rsi
	movss	(%rdi),%xmm1
	movss	(%rsi),%xmm2
	decq	%rcx
	jrcxz	inner_b_end

	# 1 left
	pslldq	$4,%xmm1
	pslldq	$4,%xmm2
	addq	$0x4,%rdi
	addq	$0x4,%rsi
	movss	(%rdi),%xmm1
	movss	(%rsi),%xmm2

inner_b_end:
	addq	$0x4,%rdi
	mulps	%xmm2,%xmm1
	addps	%xmm1,%xmm0

inner_b_0:
	# final add
	movaps	%xmm0,%xmm1
	psrldq	$8,%xmm1
	addps	%xmm1,%xmm0
	movaps	%xmm0,%xmm1
	psrldq	$4,%xmm1
	addss	%xmm1,%xmm0

	# get the "a" value and substract and divide
	movss	(%rdi),%xmm1
	subss	%xmm0,%xmm1
	divss	%xmm3,%xmm1
	movss	%xmm1,(%rdi)

	# end of alpha loop
	incq	%rbx
	cmpq	%rbx,%rdx
	jne	alpha_loop

	# end of outer loop
	incq	%rax
	cmpq	%rax,%rdx
	jne	outer_loop

	# return to calling function

	popq	%rbp
	popq	%rdi
	popq	%rsi

	popq	%r11
	popq	%r10
	popq	%r9
	popq	%r8

	popq	%rdx
	popq	%rcx
	popq	%rbx

	ret
