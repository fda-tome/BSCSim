	.file	"micro-benchmark.c"
	.text
	.section	.rodata
.LC1:
	.string	"%lu\n"
	.text
	.globl	main
	.type	main, @function
main:
.LFB6:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	subq	$80, %rsp
	movq	%fs:40, %rax
	movq	%rax, -8(%rbp)
	xorl	%eax, %eax
	movl	$0, -68(%rbp)
	movl	$400000000, %edi
	call	malloc@PLT
	movq	%rax, -56(%rbp)
	movl	$100, -60(%rbp)
	jmp	.L2
.L7:
	leaq	-48(%rbp), %rax
	movq	%rax, %rsi
	movl	$4, %edi
	call	clock_gettime@PLT
	pxor	%xmm0, %xmm0
	movss	%xmm0, -64(%rbp)
	jmp	.L3
.L4:
	call	rand@PLT
	movl	%eax, %ecx
	movl	-68(%rbp), %eax
	leal	1(%rax), %edx
	movl	%edx, -68(%rbp)
	cltq
	leaq	0(,%rax,4), %rdx
	movq	-56(%rbp), %rax
	addq	%rdx, %rax
	pxor	%xmm0, %xmm0
	cvtsi2ssl	%ecx, %xmm0
	ovss	%xmm0, (%rax)
.L3:
	movl	-68(%rbp), %eax
	cmpl	-60(%rbp), %eax
	jl	.L4
	movl	$0, -68(%rbp)
	jmp	.L5
.L6:
	movl	-68(%rbp), %eax
	leal	1(%rax), %edx
	movl	%edx, -68(%rbp)
	cltq
	leaq	0(,%rax,4), %rdx
	movq	-56(%rbp), %rax
	addq	%rdx, %rax
	movss	(%rax), %xmm0
	movss	-64(%rbp), %xmm1
	addss	%xmm1, %xmm0
	movss	%xmm0, -64(%rbp)
.L5:
	movl	-68(%rbp), %eax
	cmpl	-60(%rbp), %eax
	jl	.L6
	movl	$0, -68(%rbp)
	leaq	-32(%rbp), %rax
	movq	%rax, %rsi
	movl	$4, %edi
	call	clock_gettime@PLT
	movq	-32(%rbp), %rdx
	movq	-48(%rbp), %rax
	subq	%rax, %rdx
	imulq	$1000000, %rdx, %rsi
	movq	-24(%rbp), %rdx
	movq	-40(%rbp), %rax
	movq	%rdx, %rcx
	subq	%rax, %rcx
	movabsq	$2361183241434822607, %rdx
	movq	%rcx, %rax
	imulq	%rdx
	sarq	$7, %rdx
	movq	%rcx, %rax
	sarq	$63, %rax
	subq	%rax, %rdx
	leaq	(%rsi,%rdx), %rax
	movq	%rax, %rsi
	leaq	.LC1(%rip), %rax
	movq	%rax, %rdi
	movl	$0, %eax
	call	printf@PLT
	movl	-60(%rbp), %edx
	movl	%edx, %eax
	sall	$2, %eax
	addl	%edx, %eax
	addl	%eax, %eax
	movl	%eax, -60(%rbp)
.L2:
	cmpl	$100000000, -60(%rbp)
	jle	.L7
	nop
	movq	-8(%rbp), %rax
	subq	%fs:40, %rax
	je	.L8
	call	__stack_chk_fail@PLT
.L8:
	leave
	.cfi_def_cfa 7, 8
	ret
	.cfi_endproc
.LFE6:
	.size	main, .-main
	.ident	"GCC: (GNU) 12.1.0"
	.section	.note.GNU-stack,"",@progbits
