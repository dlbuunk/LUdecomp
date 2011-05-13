	# LU decomposition of a matrix by D.L.Buunk
	# for x86-64 processors with SSE2

	.text
	.global LUdecomp2
LUdecomp2:
	# save registers
#	pushq	%rax
	pushq	%rbx
	pushq	%rcx
	pushq	%rdx

	pushq	%r8
	pushq	%r9
	pushq	%r10
	pushq	%r11

	pushq	%r12
	pushq	%r13
	pushq	%r14
	pushq	%r15

	pushq	%rsi
	pushq	%rdi
	pushq	%rbp
	# make temp space below stack
	subq	$0x10,%rsp 	# 0x18 bytes, as to preserve alignment.

	# set up variables in their proper registers
	movq	%rdx,%rbp	# * permutations
	movq	%rcx,%rdx	# dimensions

	# set up permutation array (zero-based)
	xorq	%rax,%rax
	movq	%rdx,%rcx	# dimensions
	movq	%rbp,%rbx	# * permutations
setup_perm:
	movq	%rax,(%rbx)
	addq	$0x8,%rbx
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
	shlq	$2,%rax		# * matrix bT
	xorq	%rbx,%rbx	# size of a row of bT
setup_b:
	addq	$0x10,%rbx

	movq	%rax,(%rsi)
	addq	$0x8,%rsi
	addq	%rbx,%rax
	decq	%rdx
	jrcxz	setup_end

	movq	%rax,(%rsi)
	addq	$0x8,%rsi
	addq	%rbx,%rax
	decq	%rdx
	jrcxz	setup_end

	movq	%rax,(%rsi)
	addq	$0x8,%rsi
	addq	%rbx,%rax
	decq	%rdx
	jrcxz	setup_end

	movq	%rax,(%rsi)
	addq	$0x8,%rsi
	addq	%rbx,%rax
	decq	%rdx
	loop	setup_b
setup_end:
	# okay, setup is done, now we can setup the outer loop




	movq	%r8,%rax

	addq	$0x10,%rsp
	popq	%rbp
	popq	%rdi
	popq	%rsi

	popq	%r15
	popq	%r14
	popq	%r13
	popq	%r12

	popq	%r11
	popq	%r10
	popq	%r9
	popq	%r8

	popq	%rdx
	popq	%rcx
	popq	%rbx
#	popq	%rax

	ret
