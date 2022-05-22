# DO NOT REMOVE THIS LINE WHEN SUBMITTING
.import input.S

.text	# code starts here

.globl main
main:
	la a2 str1				# a2 = &str1
	la a3 str2 				# a3 = &str2
	lw a4 len1 				# a4 = len1
	lw a5 len2 				# a5 = len2
	jal func				# jump to func
	li a0 1					# a0 = 1
	ecall					# print a1
	li a0 10				# a0 = 10
	ecall					# exit program

func:
	li a1 0 				# res = a1 = 0
	li t0 0 				# i = t0 = 0
	addi sp sp -4			# save space for ra in stack
	sw ra 0(sp)				# store ra in stack
	jal ra loop1 			# save return address in ra
	lw ra 0(sp)				# load ra back from stack
	addi sp sp 4			# free space from stack
    jr ra 					# jump back to ra

loop1:
	bge t0 a4 done_loop1 	# if i >= len1, jump to done_loop1
	li t1 0					# j = t1 = 0
	addi sp sp -4			# save space for ra in stack
	sw ra 0(sp)				# store ra in stack
	jal ra loop2 			# save return address in ra
	lw ra 0(sp)				# load ra back from stack
	addi sp sp 4			# free space from stack
	addi t0 t0 1			# i++
	j loop1					# jump to loop1

done_loop1:
	jr ra					# jump back to ra

loop2:
	bge t1 a5 done_loop2	# if j >= len2, jump to done_loop2
	li t2 0					# k = t2 = 0
	addi sp sp -4			# save space for ra in stack
	sw ra 0(sp)				# store ra in stack
	jal ra loop3 			# save return address in ra
	lw ra 0(sp)				# load ra back from stack
	addi sp sp 4			# free space from stack
	addi t1 t1 1			# j++				
	mv a0 t2				# set a0 = t2 as parameter of get_max
	addi sp sp -4			# save space for ra in stack
	sw ra 0(sp)				# store ra in stack
	jal ra get_max			# set the larger one of k(t2 = a0) and res(a1) to a1 -> jump to get_max and set ra as the return address
	lw ra 0(sp)				# load ra back from stack
	addi sp sp 4			# free space from stack
	j loop2					# jump to loop2

done_loop2:
	jr ra					#jump back to ra

loop3:
	add t3 t0 t2 			# i+k = t3
	bge t3 a4 done_loop3	# if i+k >= len1, jump to done_loop3
	add t4 t1 t2 			# j+k = t4
	bge t4 a5 done_loop3	# if j+k >= len2, jump to done_loop3
	slli t3 t3 2 			# t3 *= 4
	slli t4 t4 2 			# t4 *= 4
	add t5 a2 t3 			# get the address of str1[i+k]
	add t6 a3 t4 			# get the address of str2[j+k]
	lw t5 0(t5) 			# t5 = str1[i+k]
	lw t6 0(t6) 			# t6 = str2[j+k]
	bne t5 t6 done_loop3 	# If str1[i+k] != str2[j+k] jump to done_loop3
	addi t2 t2 1 			# k++
	j loop3 				# jump to loop3


done_loop3:
	jr ra 					# jump back to ra

get_max:
	bge a1 a0 end			# compare the value of two parameters
	mv a1 a0				# if a1 <= a0, let a1 = a0

end:
	jr ra					# if a1 > a0, just return to ra