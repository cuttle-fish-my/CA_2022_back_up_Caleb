#   Project 1.2 RVC instructions to RISC-V instructions in RISC-V

#   You can suppose you will get the number of machine codes you should translate from 
#   "lines_of_rvc_codes" and get the 16 or 32 bits-machine codes from "rvc_codes". 
#   Both of them are stored in file "input.S"

#   When we test your project, we will change the input files content.

# DO NOT REMOVE THIS LINE WHEN SUBMITTING
.import input.S

# Your code to convert the RVC codes to RVI codes.

main:
	la a1 lines_of_rvc_codes 
	lw a1 0(a1) 		# get number of lines
	addi a1 a1 3		# prohibit overflow when calculate address
	li a2 4 		# set size of each number
	li a0 0x3cc 		# calloc
	li a6 2 		# calloc
	ecall 			# calloc
	mv s0 a0 		# save output pointer to s0
	# calloc a1*a2 bytes space for output 
	# and return pointer to a0 
	# and set all entries 0

	# allocate space for address
	li a0 0x3cc		# calloc
	ecall			# calloc
	mv s1 a0 		# save address pointer to s1

	# allocate space for bit_number
	li a0 0x3cc		# calloc
	ecall			# calloc
	mv s2 a0		# save bit_number pointer to s2
	addi a1 a1 -3		# recover a1


read_beg:
	# read instruction
	mv t2 a1 		# boundary
	li t0 0			# iterator
	li t5 0x0003 		# mask
	la s3 rvc_codes		# get address of input_instruction
read:
	bge t0 t2 read_end	# if t0 >= number of lines --> end
	lhu t3 0(s3)		# load 16-bit instruction
	and t4 t3 t5		# get last 2 bits of t3
	bne t4 t5 updReadArg_16	
	# If the last 2 bits of t3 is not 11 --> 16 bit instruction
	# If the last 2 bits of t3 is 11 --> 32 bit instruction
updReadArg_32:
	lw t3 0(s3)		# reload t3 with "load word"
	li t6 2			# t6 is the number we need to store
	li s4 0			# s4 is the iterator of bit_number array
	slli t1 t0 2		# update the offset
	add s4 s2 t1		# update the position of bit_number pointer
	sw t6 0(s4)		# store t6 in s2

	li s4 0			# s4 is the iterator of address pointer
	add s4 s1 t1		# update position of address pointer
	lw s5 0(s4)		# load address[t0] to s5
	addi s5 s5 2		# s5 = s5 + 2
	sw s5 4(s4)		# address[t0+1] = s5

	li s4 0			# s4 is the iterator of output
	add s4 s0 t1		# update position of output pointer
	sw t3 0(s4)		# output[t0] = t3

	addi t0 t0 1		# update iterator
	addi s3 s3 4		# update s3 with 32-bit

	j read			# next iteration
updReadArg_16:
	li t6 1			# t6 is the nubmer we need to store
	li s4 0			# s4 is the iterator of bit_number array
	slli t1 t0 2		# update the offset
	add s4 s2 t1		# update position of bit_number pointer
	sw t6 0(s4)		# store t6 in s2

	li s4 0			# s4 is the iterator of address pointer
	add s4 s1 t1		# update position of address pointer
	lw s5 0(s4)		# load address[t0] to s5
	addi s5 s5 1		# s5 = s5 + 1
	sw s5 4(s4)		# address[t0+1] = s5

	addi sp sp 16		# save space for t0 t1 s3
	sw t0 0(sp)		# save t0
	sw t1 4(sp)		# save t1
	sw s3 8(sp)		# save s3
	sw t2 12(sp)		# save t2

	mv a2 t3		# add parameter of check_type  
	mv a3 t0		# add parameter of check_type 
	jal ra expansion	# expand 16-bit instruction 

	lw t0 0(sp)		# load back t0
	lw t1 4(sp)		# load back t1
	lw s3 8(sp)		# load back s3
	lw t2 12(sp)		# load back t2
	addi sp sp -16		# free allocated space

	addi t0 t0 1		# update iterator
	addi s3 s3 2		# update s3 with 16-bit

	j read			# next iteration
read_end:
	j updOffset_beg		# read finished

updOffset_beg:
	li t0 0			# iterator
	la t2 lines_of_rvc_codes# load address of lines
	lw t2 0(t2)		# load lines as boundary
updOffset:
	bge t0 t2 updOffset_end	# If t0 >= number of lines --> updOffset_end
	
	li s4 0			# s4 is the iterator of s0(output)
	slli t1 t0 2		# update the offset
	add s4 s0 t1		# update the position of output address
	lw s3 0(s4)		# load output[t0] to s3

	addi sp sp 4		# save space for ra
	sw ra 0(sp)		# store ra
	mv a2 t0		# a2 is the parameter of updOffset_main --> the index of lines
	jal ra updOffset_main	# jump to updOffset_main
	lw ra 0(sp)		# load back ra
	addi sp sp -4		# free allocated space

	sw s3 0(s4)		# store the modified value back to output[t0]
	addi t0 t0 1		# update iterator

	j updOffset		# jump to updOffset

updOffset_main:
	li t3 0x0000007f	# mask of 7-bit-opcode --> 0000 0000 0000 0000 0000 0111 1111 
	and s5 s3 t3		# s5 is 7-bit-opcode 
	
	li t4 0x0000006f	# 0000 0000 0000 0000 0000 0000 0110 1111
	beq s5 t4 updOffset_jal	# --> jal

	li t4 0x00000063	# 0000 0000 0000 0000 0000 0000 0110 0011
	beq s5 t4 updOffset_SB	# --> beq or bne or blt or bge or bltu or bgeu

	ret



updOffset_end:
	j final

updOffset_SB:
	li t3 0x00000f00	# mask of imm[4:1] --> 0000 0000 0000 0000 0000 1111 0000 0000
	and s5 s3 t3		# s5 is imm[4:1] << 8
	srli s5 s5 8		# s5 is imm[4:1]

	li t3 0x7e000000	# mask of imm[10:5] --> 0111 1110 0000 0000 0000 0000 0000 0000
	and s6 s3 t3		# s6 is imm[10:5] << 21
	srli s6 s6 21		# s6 is imm[10:5]
	add s5 s5 s6		# s5 is imm[10:1]

	li t3 0x00000080	# mask of imm[7] --> 0000 0000 0000 0000 0000 0000 1000 0000
	and s6 s3 t3		# s6 is imm[7] >> 3
	slli s6 s6 3		# s6 is imm[7]
	add s5 s5 s6		# s5 is imm[11:1]

	li t3 0x80000000	# mask of imm[12] --> 1000 0000 0000 0000 0000 0000 0000 0000
	and s6 s3 t3		# s6 is imm[12] << 20
	bne s6 x0 SB_imm_neg	# If imm[12] != 0 --> imm is negative
	beq s6 x0 SB_imm_pos	# If imm[12] == 0 --> imm is positive
SB_imm_pos:
	mv t3 a2		# iterator of address(s1)

SB_imm_pos_loops_beg:
	ble s5 x0 SB_imm_pos_loops_end
	li t6 0			# t6 is iterator_pointer of address(s1)
	slli t4 t3 2		# update offset

	add t6 s2 t4		# update position of t6
	lw t6 0(t6)		# t6 = address[t3]

	sub s5 s5 t6		# update imm (or offset)
	addi t3 t3 1		# update iterator
	mv a0 t3		# set return value a0 = t3

	j SB_imm_pos_loops_beg	# next iteration

SB_imm_pos_loops_end:
	sub t3 a0 a2 		# get difference of line index
	mv t5 t3
	add t5 t5 t5
	
	li t3 0x01fff07f	# mask of imm[12:1] --> 0000 0001 1111 1111 1111 0000 0111 1111
	and s3 s3 t3		# the imm of s3 is setted to all 0
	
	li t3 0x0000000f	# mask of imm[4:1] --> 0000 0000 0000 0000 0000 0000 0000 1111
	and t6 t5 t3		# t6 is imm[4:1] >> 8
	slli t6 t6 8		# t6 is imm[4:1]
	add s3 s3 t6		# finish imm[4:1]
	
	li t3 0x000003f0	# mask of imm[10:5] --> 0000 0000 0000 0000 0000 0011 1111 0000
	and t6 t5 t3		# t6 is imm[10:5] >> 21
	slli t6 t6 21		# t6 is imm[10:5]
	add s3 s3 t6		# finish imm[10:5]

	li t3 0x00000400	# mask of imm[11] --> 0000 0000 0000 0000 0000 0100 0000 0000
	and t6 t5 t3		# t6 is imm[11] << 3
	srli t6 t6 3		# t6 is imm[11]
	add s3 s3 t6		# finish imm[11]

	j updOffset_SB_end

SB_imm_neg:
	srli s6 s6 20		# t6 is imm[12]
	add s5 s5 s6		# s5 is imm[12:1]
	slli s5 s5 20		# make s5 two-complemented 
	srai s5 s5 20		# make s5 two-complemented

	mv t3 a2		# iterator of address(s1)
	addi t3 t3 -1		# we do not count the line where we are
	
SB_imm_neg_loops_beg:
	bge s5 x0 SB_imm_neg_loops_end 
	li t6 0			# t6 is iterator_pointer of address(s1)
	slli t4 t3 2		# update offset

	add t6 s2 t4		# update position of t6
	lw t6 0(t6)		# t6 = address[t3]

	add s5 s5 t6		# update imm (or offset)
	mv a0 t3		# set return value a0 = t3
	addi t3 t3 -1		# update iterator

	j SB_imm_neg_loops_beg	# next iteration

SB_imm_neg_loops_end:
	sub t3 a0 a2 		# get difference of line index
	mv t5 t3
	add t5 t5 t5
	
	li t3 0x01fff07f	# mask of imm[12:1] --> 0000 0001 1111 1111 1111 0000 0111 1111
	and s3 s3 t3		# the imm of s3 is setted to all 0
	
	li t3 0x0000000f	# mask of imm[4:1] --> 0000 0000 0000 0000 0000 0000 0000 1111
	and t6 t5 t3		# t6 is imm[4:1] >> 8
	slli t6 t6 8		# t6 is imm[4:1]
	add s3 s3 t6		# finish imm[4:1]
	
	li t3 0x000003f0	# mask of imm[10:5] --> 0000 0000 0000 0000 0000 0011 1111 0000
	and t6 t5 t3		# t6 is imm[10:5] >> 21
	slli t6 t6 21		# t6 is imm[10:5]
	add s3 s3 t6		# finish imm[10:5]

	li t3 0x00000400	# mask of imm[11] --> 0000 0000 0000 0000 0000 0100 0000 0000
	and t6 t5 t3		# t6 is imm[11] << 3
	srli t6 t6 3		# t6 is imm[11]
	add s3 s3 t6		# finish imm[11]

	li t3 0x80000000	# t3 is imm[12]
	add s3 s3 t3		# finish imm[12]

	j updOffset_SB_end

updOffset_SB_end:
	ret

updOffset_jal:
	# update imm of jal (a2 is the index of lines and s1 is the address)
	# your code here
	li t3 0x7fe00000	# mask of imm[10:1] --> 0111 1111 1110 0000 0000 0000 0000 0000
	and t5 s3 t3		# t5 is imm[10:1] << 21
	srli t5 t5 21		# t5 is imm[10:1]

	li t3 0x00100000	# mask of imm[11] --> 0000 0000 0001 0000 0000 0000 0000 0000 	
	and t6 s3 t3		# t6 is imm[11] << 10
	srli t6 t6 10		# t6 is imm[11]
	add t5 t5 t6		# t5 is imm[11:1]

	li t3 0x000ff000	# mask of imm[19:12] --> 0000 0000 0000 1111 1111 0000 0000 0000
	and t6 s3 t3		# t6 is imm[19:12] << 1
	srli t6 t6 1		# t6 is imm[19:12]
	add t5 t5 t6		# t5 is imm[19:20]

	li t3 0x80000000	# mask of imm[20] --> 1000 0000 0000 0000 0000 0000 0000 0000
	and t6 s3 t3		# t6 is imm[20] << 12
	bne t6 x0 jal_imm_neg	# If imm[20] << 12 == 0 --> imm is positive
	beq t6 x0 jal_imm_pos	# If imm[20] << 12 != 0 --> imm is negative
jal_imm_pos:
	mv t3 a2		# iterator of address(s1)
jal_imm_pos_loops_beg:
	ble t5 x0 jal_imm_pos_loops_end
	li t6 0			# t6 is iterator_pointer of address(s1)
	slli t4 t3 2 		# update offset

	add t6 s2 t4		# update position of t6
	lw t6 0(t6)		# t6 = address[t3]

	sub t5 t5 t6		# update imm (or offset) 
	addi t3 t3 1		# update iterator
	mv a0 t3		# set return value a0 = t3

	j jal_imm_pos_loops_beg

jal_imm_pos_loops_end:
	sub t3 a0 a2		# get difference of line index
	mv t5 t3
	add t5 t5 t5

	li t3 0x00000fff	# mask of imm[20:1]
	and s3 s3 t3		# the imm of s3 is set to all 0

	li t3 0x000003ff	# mask of imm[10:1] --> 0000 0000 0000 0000 0000 0011 1111 1111
	and t6 t5 t3		# t6 is imm[10:1] >> 21
	slli t6 t6 21		# t6 is imm[10:1]
	add s3 s3 t6		# finish imm[10:1]

	li t3 0x00000400	# mask of imm[11] --> 0000 0000 0000 0000 0000 0100 0000 0000
	and t6 t5 t3		# t6 is imm[11] >> 10
	slli t6 t6 10		# t6 is imm[11]
	add s3 s3 t6		# finish imm[11]

	li t3 0x0007f800 	# mask of imm[19:12] --> 0000 0000 0000 0111 1111 1000 0000 0000
	and t6 t5 t3		# t6 is imm[19:12] >> 1
	slli t6 t6 1		# t6 is imm[19:12]
	add s3 s3 t6		# finish imm[19:12]

	# imm[20] = 0 --> no need to update

	j updOffset_jal_end


jal_imm_neg:
	srai t6 t6 12 		# t6 is imm[20]
	add t5 t5 t6		# t5 is imm[20:1]

	mv t3 a2		# iterator of address(s1)	

	addi t3 t3 -1		# we do not count the line where we just stay

jal_imm_neg_loops_beg:
	bge t5 x0 jal_imm_neg_loops_end 
	li t6 0			# t6 is iterator_pointer of address(s1)
	slli t4 t3 2		# update offset

	add t6 s2 t4		# update position of t6
	lw t6 0(t6)		# t6 = address[t3]

	add t5 t5 t6		# update imm (or offset)
	mv a0 t3		# set return value a0 = t3
	addi t3 t3 -1		# update iterator

	j jal_imm_neg_loops_beg	# next iteration
	
jal_imm_neg_loops_end:
	sub t3 a0 a2		# get difference of line index
	mv t5 t3
	add t5 t5 t5

	li t3 0x00000fff	# mask of imm[20:1]
	and s3 s3 t3		# the imm of s3 is setted to all 0

	li t3 0x000003ff	# mask of imm[10:1] --> 0000 0000 0000 0000 0000 0011 1111 1111
	and t6 t5 t3		# t6 is imm[10:1] >> 21
	slli t6 t6 21		# t6 is imm[10:1]
	add s3 s3 t6		# finish imm[10:1]

	li t3 0x00000400	# mask of imm[11] --> 0000 0000 0000 0000 0000 0100 0000 0000
	and t6 t5 t3		# t6 is imm[11] >> 10
	slli t6 t6 10		# t6 is imm[11]
	add s3 s3 t6		# finish imm[11]

	li t3 0x0007f800 	# mask of imm[19:12] --> 0000 0000 0000 0111 1111 1000 0000 0000
	and t6 t5 t3		# t6 is imm[19:12] >> 1
	slli t6 t6 1		# t6 is imm[19:12]
	add s3 s3 t6		# finish imm[19:12]

	li t3 0x80000000	# t3 imm[20] --> 1000 0000 0000 0000 0000 0000 0000 0000
	add s3 s3 t3		# finish imm[20]

	j updOffset_jal_end


updOffset_jal_end:
	ret

print_beg:
	li s7 0			# iterator
	la s8 lines_of_rvc_codes# get address of lines
	lw s8 0(s8)		# boundary
	mv s6 s0 		# s6 is the to-be tested array
	# test s0 --> lw s6 0(sp)
	# test s1 --> lw s6 4(sp)
	# test s2 --> lw s6 8(sp)
print:
	bge s7 s8 print_end	# if s7 >= lines --> end
	lw a2 0(s6)		# load value of entry s6 pointing to
	li t0 31		# iterator
	li t1 0			# boundary
	li t2 1			# mask
	slli t2 t2 31		# 0x80000000

print_bit:
	blt t0 t1 print_bit_end
	and t3 a2 t2		
	srl t3 t3 t0
	addi t0 t0 -1
	srli t2 t2 1 
	li a0 1
	mv a1 t3
	ecall
	j print_bit

print_bit_end:
	li a1 10		# set '\n' between entries
	li a0 11		# print character
	ecall			# print character
	addi s6 s6 4		# let s6 point to the next entry
	addi s7 s7 1		# update iterator
	j print			# loops
print_end:
	ret

final:
	addi sp sp 4		# save space for ra
	sw ra 0(sp)		# save ra in 0(sp)
	jal ra print_beg	# jump to print
	lw ra 0(sp)		# load back ra from 0(sp)
	addi sp sp -4		# free allocated space of sp
end:
	# Exit with error code 0
	add a1, x0, x0
	addi a0, x0, 17
	ecall


expansion:
	addi sp sp 4		
	sw ra 0(sp)
	jal ra check_type 	# check the type of instruction and return the corresponding code
	lw ra 0(sp)
	addi sp sp -4
	# Here is the expansion code
	# Here is the expansion code
	li s9 1    	        #type=1
	beq a0 s9 expand_c_add  #expand_c_add
	li s9 2			#type=2
	beq a0 s9 expand_c_mv 	#expand_c_mv
	li s9 3			#type=3
	beq a0 s9 expand_c_jr 	#expand_c_jr
	li s9 4			#type=4
	beq a0 s9 expand_c_jalr	#expand_c_jalr
	li s9 5			#type=5
	beq a0 s9 expand_c_li   #expand_c_li
	li s9 6			#type=6
	beq a0 s9 expand_c_lui 	#expand_c_lui
	li s9 7			#type=7
	beq a0 s9 expand_c_addi	#expand_c_addi
	li s9 8			#type=8
	beq a0 s9 expand_c_slli	#expand_c_slli
	li s9 9			#type=9
	beq a0 s9 expand_c_lw 	#expand_c_lw
	li s9 10		#type=10
	beq a0 s9 expand_c_sw	#expand_c_sw
	li s9 11		#type=11
	beq a0 s9 expand_c_and 	#expand_c_and
	li s9 12		#type=12
	beq a0 s9 expand_c_or	#expand_c_or
	li s9 13 		#type=13
	beq a0 s9 expand_c_xor	#expand_c_xor
	li s9 14 		#type=14
	beq a0 s9 expand_c_sub	#expand_c_sub
	li s9 15		#type=15
	beq a0 s9 expand_c_beqz_tmp	#expand_c_beqz
	li s9 16		#type=16
	beq a0 s9 expand_c_bnez_tmp	#expand_c_bneq
	li s9 17 		#type=17
	beq a0 s9 expand_c_srli_tmp #expand_c_srli
	li s9 18		#type=18
	beq a0 s9 expand_c_srai_tmp	#expand_c_srai
	li s9 19		#type=19
	beq a0 s9 expand_c_andi_tmp 	#expand_c_andi
	li s9 20		#type=20
	beq a0 s9 expand_c_j_tmp	#expand_c_j
	li s9 21		#type=21
	beq a0 s9 expand_c_jal_tmp 	#expand_c_jal
	# If a0 == 1 to 21
	# doing the corresponding expansion
	# the format of label is expand_{type}

	# Here is the expansion code
	# Here is the expansion code

	Exit:
	ret

check_type:
	li t0 0x0003 		# mask of opcode --> 0000 0000 0000 0011
	and s3 a2 t0		# s3 is opcode

	li t1 0x0002		# 0000 0000 0000 0010
	beq s3 t1 check_CR_CI 	# --> CR or CI type

	li t1 0x0001		# 0000 0000 0000 0001
	beq s3 t1 check_CI_CA_CB_CJ	# --> CI or CA or CB or Cj type

	li t1 0x0000		# 0000 0000 0000 0000
	beq s3 t1 check_CL_CS	# --> CL or CS type

	ret


check_CR_CI:
	li t0 0xf000 		# mask of funct4 --> 1111 0000 0000 0000
	and s3 a2 t0		# s3 is funct4

	li t1 0x9000		# 1001 0000 0000 0000
	beq s3 t1 check_add_jalr# --> add or jalr

	li t1 0x8000		# 1000 0000 0000 0000
	beq s3 t1 check_mv_jr	# --> mv or jr

	li t0 0xe000		# mask of funct3 --> 1110 0000 0000 0000
	and s3 a2 t0		# s3 is funct3

	li t1 0x0000		# 0000 0000 0000 0000
	beq s3 t1 check_slli	# --> slli

	ret

check_CI_CA_CB_CJ:
	li t0 0xe000		# mask of funct3 --> 1110 0000 0000 0000
	and s3 a2 t0		# s3 is funct3

	li t1 0x4000		# 0100 0000 0000 0000
	beq s3 t1 check_li	# --> li

	li t1 0x6000		# 0110 0000 0000 0000
	beq s3 t1 check_lui	# --> lui

	li t1 0x0000		# 0000 0000 0000 0000
	beq s3 t1 check_addi# --> addi

	li t1 0xc000		# 1100 0000 0000 0000
	beq s3 t1 check_beqz# --> beqz

	li t1 0xe000		# 1110 0000 0000 0000
	beq s3 t1 check_bnez# --> bnez

	li t1 0x8000		# 1000 0000 0000 0000
	beq s3 t1 check_srli_srai_andi	# --> srli or srai or andi

	li t1 0xa000		# 1010 0000 0000 0000
	beq s3 t1 check_j	# --> j

	li t1 0x2000		# 0010 0000 0000 0000
	beq s3 t1 check_jal	# --> jal
check_CA:
	li t0 0x0060		# mask of funct2 --> 0000 0000 0110 0000
	and s3 a2 t0		# s3 is funct2

	li t1 0x0060		# 0000 0000 0110 0000
	beq s3 t1 check_and 	# --> and

	li t1 0x0040		# 0000 0000 0100 0000
	beq s3 t1 check_or 	# --> or

	li t1 0x0020		# 0000 0000 0010 0000
	beq s3 t1 check_xor	# --> xor

	li t1 0x0000		# 0000 0000 0000 0000
	beq s3 t1 check_sub	# --> sub

	ret
check_CL_CS:
	li t0 0xe000		# mask of funct3 --> 1110 0000 0000 0000
	and s3 a2 t0		# s3 is funct3

	li t1 0x4000		# 0100 0000 0000 0000
	beq s3 t1 check_lw	# --> lw

	li t1 0xc000		# 1100 0000 0000 0000
	beq s3 t1 check_sw	# --> sw

	ret


check_add_jalr:
	li t0 0x007c		# mask of rs2
	and s3 a2 t0		# s3 is rs2

	beqz s3 check_jalr	# --> jalr
	bnez s3 check_add	# --> add

	ret
check_mv_jr:
	li t0 0x007c		# mask of rs2
	and s3 a2 t0		# s3 is rs2

	beqz s3 check_jr	# --> jr
	bnez s3 check_mv	# --> mv

	ret

check_srli_srai_andi:
	li t0 0x0c00		# mask of funct2 --> 0000 1100 0000 0000
	and s3 a2 t0		# s3 is funct2

	li t1 0x0000		# 0000 0000 0000 0000
	beq s3 t1 check_srli	# --> srli

	li t1 0x0400		# 0000 0100 0000 0000
	beq s3 t1 check_srai	# --> srai

	li t1 0x0800		# 0000 1000 0000 0000
	beq s3 t1 check_andi	# --> andi

	j check_CA

check_add:
	li a0 1
	ret
check_mv:
	li a0 2
	ret
check_jr:
	li a0 3
	ret
check_jalr:
	li a0 4
	ret
check_li:
	li a0 5
	ret
check_lui:
	li a0 6
	ret
check_addi:
	li a0 7
	ret
check_slli:
	li a0 8
	ret
check_lw:
	li a0 9
	ret
check_sw:
	li a0 10
	ret
check_and:
	li a0 11
	ret
check_or:
	li a0 12
	ret
check_xor:
	li a0 13
	ret
check_sub:
	li a0 14
	ret
check_beqz:
	li a0 15
	ret
check_bnez:
	li a0 16
	ret
check_srli:
	li a0 17
	ret
check_srai:
	li a0 18
	ret
check_andi:
	li a0 19
	ret
check_j:
	li a0 20
	ret
check_jal:
	li a0 21
	ret
expand_c_beqz_tmp:
	j expand_c_beqz
expand_c_bnez_tmp:
	j expand_c_bnez
expand_c_srli_tmp:
	j expand_c_srli
expand_c_srai_tmp:
	j expand_c_srai
expand_c_andi_tmp:
	j expand_c_andi
expand_c_j_tmp:
	j expand_c_j
expand_c_jal_tmp:
	j expand_c_jal
expand_c_add:
	#a2 16 bits word
	#a3 iterator
	 
	srli a2 a2 2		# erase OPcode
	li t0 0x1F		# RS2 filter 11111
        and t1 a2 t0		# t1 rs2  
	li t0 0x3E0		# RD filter 1111100000
	and t2 a2 t0		# t2 rd|00000
	srli t2 t2 5            # t2 rd
	li t3 0x33 		# t3 OPcode
	slli t2 t2 7            # rd|0000000 
	add t3 t3 t2 		# rd|OPcode
	slli t2 t2 8		# rd|000|00000|0000000
	add t3 t3 t2		# rd|func3|rd|OPcode
	slli t1 t1 20		# rs2|00000|000|00000|0000000
	add t3 t3 t1		# func7|rs2|rd|func3|rd|OPcode

	slli t0 a3 2 		# update the offset
	li s4 0			# s4 is the iterator of output
	add s4 s0 t0		# update position of output pointer
	sw t3 0(s4)		# output[t0] = t3
	j Exit

expand_c_mv:
	#a2 16 bits word
	#a3 iterator
	 
	srli a2 a2 2		# erase OPcode
	li t0 0x1F		# RS2 filter 11111
        and t1 a2 t0		# t1 rs2  
	li t0 0x3E0		# RD filter 1111100000
	and t2 a2 t0		# t2 rd|00000
	srli t2 t2 5            # t2 rd
	li t3 0x33 		# t3 OPcode 0110011
	slli t2 t2 7          	# rd|0000000 
	add t3 t3 t2 		# rd|OPcode
	slli t1 t1 20		# rs2|00000|000|00000|0000000
	add t3 t3 t1		# func7|rs2|x0|func3|rd|OPcode

	slli t0 a3 2 		# update the offset
	li s4 0			# s4 is the iterator of output
	add s4 s0 t0		# update position of output pointer
	sw t3 0(s4)		# output[t0] = t3
	j Exit

expand_c_jr:
	#a2 16 bits word
	#a3 iterator
	 
	srli a2 a2 2		# erase OPcode
	li t0 0x3E0		# RS1 filter 1111100000
	and t2 a2 t0		# t2 rs1|00000
	li t3 0x67 		# t3 OPcode 1100111
	slli t2 t2 10           # rs1|000|00000|0000000 
	add t3 t3 t2 		# 00000000000|rs1|func3|x0|OPcode

	slli t0 a3 2 		# update the offset
	li s4 0			# s4 is the iterator of output
	add s4 s0 t0		# update position of output pointer
	sw t3 0(s4)		# output[t0] = t3
	j Exit

expand_c_jalr:
	#a2 16 bits word
	#a3 iterator
	 
	srli a2 a2 2		# erase OPcode
	li t0 0x3E0		# RS1 filter 1111100000
	and t2 a2 t0		# t2 rs1|00000
	li t3 0xe7 		# t3 x1|OPcode 00001|1100111
	slli t2 t2 10           # rs1|000|00000|0000000 
	add t3 t3 t2 		# 00000000000|rs1|func3|x1|OPcode	

	slli t0 a3 2 		# update the offset
	li s4 0			# s4 is the iterator of output
	add s4 s0 t0		# update position of output pointer
	sw t3 0(s4)		# output[t0] = t3
	j Exit



expand_c_li:
	#a2 16 bits word
	#a3 iterator
	 
	srli a2 a2 2		# erase OPcode
	li t0 0x1F              # imm[4:0] filter 11111
	and t1 a2 t0		# t1 imm[4:0]

	li t0 0x3E0		# RD filter 1111100000
	and t2 a2 t0		# t2 rd|00000

	li t3 0x13 		# t3 OPcode 0010011
	slli t2 t2 2            # rd|0000000 
	add t3 t3 t2 		# x0|func3|rd|OPcode
	
        slli t1 t1 20           # imm[4:0]|00000|000|00000|0000000
	add t3 t3 t1		# imm[4:0]|x0|func3|rd|opcode

	slli a2 a2 2
	li t0 0x1000            # imm[5] filter 1000000000000
	and t1 a2 t0 		# t1 imm[5]|12X0
	slli t1 t1 19		# imm[5]|31X0
	srai t1 t1 6
	add t3 t3 t1 		# imm[5]...imm[4:0]|x0|func3|rd|OPcode

	slli t0 a3 2 		# update the offset
	li s4 0			# s4 is the iterator of output
	add s4 s0 t0		# update position of output pointer
	sw t3 0(s4)		# output[t0] = t3
	j Exit

expand_c_lui:
	#a2 16 bits word
	#a3 iterator
	 
	srli a2 a2 2		# erase OPcode
	li t0 0x1F              # nzimm[16:12] filter 11111
	and t1 t1 t0		# t1 nzimm[16:12]

	and t1 a2 t0		# t1 nzimm[16:12]
	li t0 0x3E0		# RD filter 1111100000
	and t2 a2 t0		# t2 rd|00000
	srli t2 t2 5		# t2 rd

	li t3 0x37 		# t3 OPcode 0110111
	slli t2 t2 7            # rd|0000000 
	add t3 t3 t2 		# rd|OPcode

        slli t1 t1 12           # nzimm[16:12]|00000|0000000
	add t3 t3 t1		# nzimm[16:12]|rd|opcode

	li t0 0x400          	# nzimm[17] filter 1|00000|00000
	and t1 a2 t0 		# t1 nzimm[17]|12X0
	srli t1 t1 10

	slli t1 t1 31		# nzimm[17]...
	srai t1 t1 14
	add t3 t3 t1 		# nzimm[17]...nzimm[16:12]|rd|OPcode

	slli t0 a3 2 		# update the offset
	li s4 0			# s4 is the iterator of output
	add s4 s0 t0		# update position of output pointer
	sw t3 0(s4)		# output[t0] = t3
	j Exit

expand_c_addi:
	#a2 16 bits word
	#a3 iterator
	 
	srli a2 a2 2		# erase OPcode
	li t0 0x1F              # nzimm[4:0] filter 11111
	and t1 a2 t0		# t1 nzimm[4:0]

	li t0 0x3E0		# RD filter 1111100000
	and t2 a2 t0		# t2 rd|00000
	srli t2 t2 5		# t2 rd

	li t3 0x13 		# t3 OPcode 0010011
	slli t2 t2 7            # rd|0000000 
	add t3 t3 t2 		# rd|OPcode

	slli t2 t2 8            # rd|000|00000|0000000 
	add t3 t3 t2 		# rd|func3|rd|OPcode

        slli t1 t1 20           # nzimm[4:0]|00000|000|00000|0000000
	add t3 t3 t1		# nzimm[4:0]|rd|func3|rd|opcode

	li t0 0x400            # nzimm[5] filter 10000000000
	and t1 a2 t0 		# t1 nzimm[5]|10X0
	srli t1 t1 10            

	slli t1 t1 31		# nzimm[5]...
	srai t1 t1 6
	add t3 t3 t1 		# nzimm[5]...imm[4:0]|rd|func3|rd|OPcode

	slli t0 a3 2 		# update the offset
	li s4 0			# s4 is the iterator of output
	add s4 s0 t0		# update position of output pointer
	sw t3 0(s4)		# output[t0] = t3
	j Exit

expand_c_slli:
	#a2 16 bits word
	#a3 iterator
	 
	srli a2 a2 2		# erase OPcode
	li t0 0x1F              # shamt[4:0] filter 11111
	and t1 a2 t0		# t1 shamt[4:0]

	li t0 0x3E0		# RD filter 1111100000
	and t2 a2 t0		# t2 rd|00000
	srli t2 t2 5		# t2 rd

	li t3 0x13 		# t3 OPcode 0010011
	slli t2 t2 7            # rd|0000000 
	add t3 t3 t2 		# rd|OPcode

	li t0 0x1000		# func3 001|00000|0000000
	add t3 t3 t0 		# fun3|rd|OPcode

	slli t2 t2 8            # rd|000|00000|0000000 
	add t3 t3 t2 		# rd|func3|rd|OPcode

        slli t1 t1 20           # shamt[4:0]|00000|000|00000|0000000
	add t3 t3 t1		# func7|shamt[4:0]|rd|func3|rd|opcode
	
	slli t0 a3 2 		# update the offset
	li s4 0			# s4 is the iterator of output
	add s4 s0 t0		# update position of output pointer
	sw t3 0(s4)		# output[t0] = t3
	j Exit

expand_c_lw:
	li t0 0x1C		# RD' filter 11100
	and t1 a2 t0 		# t1 RD'|00
	li t0 0x20 		# ex-RD 01000|00
	add t1 t1 t0		# t1 RD|00

	li t0 0x380		# RS' filter 1110000000
	and t2 a2 t0 		# t2 RS'|00|000|00
	li t0 0x400		# ex-RS 01000|00|000|00
	add t2 t2 t0            # RS|00|000|00

	li t3 0x3		# opcode 0000011
	slli t1 t1 5		# t1 rd|0000000
	add t3 t3 t1 		# rd|opcode

	li t0 0x2000            # 010|00000|0000000
	add t3 t3 t0 		# func3|rd|opcode

	slli t2 t2 8            # rs|000|00000|0000000
	add t3 t3 t2 		# rs|func3|rd|opcode

	li t0 0x40		# offset[2] filter 1000000
	and t1 a2 t0 		# t1 offset[2]0|000|00
	slli t1 t1 16           # offset[2]00|00000|000|00000|0000000
	add t3 t3 t1 		# offset[2]00|rs|func3|rd|opcode

	li t0 0x1C00		# offset[5:3] filter 111|000|00|000|00
	and t1 a2 t0 		# t1 offset[5:3]|000|00|000|00
	slli t1 t1 13		# offset[5:3]000|00000|000|00000|0000000
	add t3 t3 t1 		# offset[5:2]00|rs|func3|rd|opcode

	li t0 0x20		# offset[6] filter 100000
	and t1 a2 t0 		# t1 offset[6]|000|00
	slli t1 t1 21           # offset[6]000000|00000|000|00000|0000000
	add t3 t3 t1 		# offset[6:2]00|rs|func3|rd|opcode

	slli t0 a3 2 		# update the offset
	li s4 0			# s4 is the iterator of output
	add s4 s0 t0		# update position of output pointer
	sw t3 0(s4)		# output[t0] = t3
	j Exit

expand_c_sw:
	li t0 0x1C		# RS2' filter 11100
	and t1 a2 t0 		# t1 RS2'|00
	li t0 0x20 		# ex-RS2 01000|00
	add t1 t1 t0		# t1 RS2|00

	li t0 0x380		# RS1' filter 1110000000
	and t2	a2 t0 		# t2 RS1'|00|000|00
	
	li t0 0x400		# ex-RS1 01000|00|000|00
	add t2 t2 t0            # RS1|00|000|00

	li t3 0x23		# opcode 0100011
	slli t2 t2 8		# rs1|000|00000|0000000
	add t3 t3 t2 		# rs1|000|00000|opcode


	li t0 0x2000            # 010|00000|0000000
	add t3 t3 t0 		# rs1|func3|00000|opcode

	slli t1 t1 18           # rs2|00000|000|00000|0000000
	add t3 t3 t1 		# rs2|rs1|func3|00000|opcode

	li t0 0x40		# offset[2] filter 1000000
	and t1 a2 t0 		# t1 offset[2]0|000|00
	slli t1 t1 3            # offset[2]00|0000000
	add t3 t3 t1 		# rs2|rs1|func3|00offset[2]00|opcode

	li t0 0xC00		# offset[4:3] filter 11|000|00|000|00
	and t1 a2 t0 		# t1 offset[4:3]|000|00|000|00
	slli t1 t1 0		# offset[4:3]000|0000000
	add t3 t3 t1 		# rs2|rs1|func3|offset[4:2]00|opcode

	li t0 0x1000		# offset[5] filter 100|000|00|000|00
	and t1 a2 t0 		# t1 offset[5]00|000|00|000|00
	slli t1 t1 13           # offset[5]|00000|00000|000|00000|0000000
	add t3 t3 t1 		# offset[5]|rs2|rs1|func3|offset[4:2]00|opcode

	li t0 0x20		# offset[6] filter 100000
	and t1 a2 t0 		# t1 offset[6]|000|00
	slli t1 t1 21           # offset[6]0|00000|00000|000|00000|0000000
	add t3 t3 t1 		# offset[6:5]|rs2|rs1|func3|offset[4:2]00|opcode
	
	slli t0 a3 2 		# update the offset
	li s4 0			# s4 is the iterator of output
	add s4 s0 t0		# update position of output pointer
	sw t3 0(s4)		# output[t0] = t3
	j Exit


expand_c_and:
	li t0 0x1C		# RS2' filter 11100
	and t1 a2 t0 		# t1 RS2'|00
	li t0 0x20 		# ex-Rs2 01000|00
	add t1 t1 t0		# t1 Rs2|00

	li t0 0x380		# RD' filter 1110000000
	and t2 a2 t0 		# t2 RD'|00|000|00
	li t0 0x400		# ex-RD 01000|00|000|00
	add t2 t2 t0            # RD|00|000|00
	
	li t3 0x33 		# opcode 011011
	add t3 t3 t2 		# rd|opcode
	
	li t0 0x7000		# func3 111|00000|0000000
	add t3 t3 t0 		# func3|rd|opcode

	slli t2 t2 8  		# rd|000|00000|0000000
	add t3 t3 t2 		# rd|func3|rd|opcode

	slli t1 t1 18 		# rs2|00000|000|00000|0000000
	add t3 t3 t1 		# func7|rs2|rd|func3|rd|opcode

	slli t0 a3 2 		# update the offset
	li s4 0			# s4 is the iterator of output
	add s4 s0 t0		# update position of output pointer
	sw t3 0(s4)		# output[t0] = t3
	j Exit

expand_c_or:
	li t0 0x1C		# RS2' filter 11100
	and t1 a2 t0 		# t1 RS2'|00
	li t0 0x20 		# ex-RD 01000|00
	add t1 t1 t0		# t1 RD|00

	li t0 0x380		# RD' filter 1110000000
	and t2	a2 t0 		# t2 RD'|00|000|00
	li t0 0x400		# ex-RD 01000|00|000|00
	add t2 t2 t0            # RD|00|000|00

	li t3 0x33 		# opcode 011011
	add t3 t3 t2 		# rd|opcode
	
	li t0 0x6000		# func3 110|00000|0000000
	add t3 t3 t0 		# func3|rd|opcode

	slli t2 t2 8  		# rd|000|00000|0000000
	add t3 t3 t2 		# rd|func3|rd|opcode

	slli t1 t1 18 		# rs2|00000|000|00000|0000000
	add t3 t3 t1 		# func7|rs2|rd|func3|rd|opcode
	 
	slli t0 a3 2 		# update the offset
	li s4 0			# s4 is the iterator of output
	add s4 s0 t0		# update position of output pointer
	sw t3 0(s4)		# output[t0] = t3
	j Exit

expand_c_xor:
	li t0 0x1C		# RS2' filter 11100
	and t1 a2 t0 		# t1 RS2'|00
	li t0 0x20 		# ex-RD 01000|00
	add t1 t1 t0		# t1 RD|00

	li t0 0x380		# RD' filter 1110000000
	and t2	a2 t0 		# t2 RD'|00|000|00
	li t0 0x400		# ex-RD 01000|00|000|00
	add t2 t2 t0            # RD|00|000|00

	li t3 0x33 		# opcode 011011
	add t3 t3 t2 		# rd|opcode
	
	li t0 0x4000		# func3 100|00000|0000000
	add t3 t3 t0 		# func3|rd|opcode

	slli t2 t2 8  		# rd|000|00000|0000000
	add t3 t3 t2 		# rd|func3|rd|opcode

	slli t1 t1 18 		# rs2|00000|000|00000|0000000
	add t3 t3 t1 		# func7|rs2|rd|func3|rd|opcode
	 
	slli t0 a3 2 		# update the offset
	li s4 0			# s4 is the iterator of output
	add s4 s0 t0		# update position of output pointer
	sw t3 0(s4)		# output[t0] = t3
	j Exit

expand_c_sub:
	li t0 0x1C		# RS2' filter 11100
	and t1 a2 t0 		# t1 RS2'|00
	li t0 0x20 		# ex-RD 01000|00
	add t1 t1 t0		# t1 RD|00

	li t0 0x380		# RD' filter 1110000000
	and t2	a2 t0 		# t2 RD'|00|000|00
	li t0 0x400		# ex-RD 01000|00|000|00
	add t2 t2 t0            # RD|00|000|00

	li t3 0x33 		# opcode 011011
	add t3 t3 t2 		# rd|opcode
	
	slli t2 t2 8  		# rd|000|00000|0000000
	add t3 t3 t2 		# rd|func3|rd|opcode

	slli t1 t1 18 		# rs2|00000|000|00000|0000000
	add t3 t3 t1 		# rs2|rd|func3|rd|opcode

	li t0 0x40000000	# func7 0100000|00000|00000|000|00000|0000000
	add t3 t3 t0 		# func7|rs2|rd|func3|rd|opcode	

	slli t0 a3 2 		# update the offset
	li s4 0			# s4 is the iterator of output
	add s4 s0 t0		# update position of output pointer
	sw t3 0(s4)		# output[t0] = t3
	j Exit




expand_c_beqz:
	li t0 0x380 		# RS1' filter 111|00000|00
	and t2 a2 t0 		# t2 RS1'|00000|00
	li t0 0x400		# ex-RS1 01000|00000|00
	add t2 t2 t0            # RS1|00000|00
	 
	li t0 0x18   		# offset[2|1]filter 11000
	and t1 a2 t0  	        # offset[2|1]0|00

	li t0 0xC00		# offset[4|3]filter 11|000|00000|00
	and t3 a2 t0 		# offset[4|3]|000|00000|00
	slli t1 t1 5            # offset[2|1]0|0000000
	add t1 t1 t3 		# offset[4|1]0|0000000

	li t3 0x63		# opcode 1100011
	add t3 t3 t1 		# func3|offset[4:1]0|opcode

	slli t2 t2 8            # rs1|000|00000|0000000
	add t3 t3 t2 		# x0|rs1|func3|func3|offset[4:1]0|opcode

	li t0 0x4  		# offset[5]filter 100
	and t1 a2 t0 		# offset[5]|00
	slli t1 t1 23           # offset[5]|00000|00000|000|00000|0000000
	add t3 t3 t1 		# offset[5]|x0|rs1|func3|func3|offset[4:1]0|opcode

	li t0 0x60 		# offset[7|6]filter 1100000
	and t1 a2 t0            # offset[7|6]000|00
	slli t1 t1 21           # offset[7|6]0|00000|00000|000|00000|0000000
	add t3 t3 t1 	 	# offset[7:5]|x0|rs1|func3|func3|offset[4:1]0|opcode

	li t0 0x1000 		# offset[8]filter 100|000|00000|00
	and t1 a2 t0            # offset[8]00|000|00000|00
	srli t1 t1 5
	add t3 t3 t1
	slli t1 t1 24           # offset[8]|0X31
	srai t1 t1 3
	add t3 t3 t1 	 	# offset[8]...offset[7:5]|x0|rs1|func3|func3|offset[4:1]0|opcode


	slli t0 a3 2 		# update the offset
	li s4 0			# s4 is the iterator of output
	add s4 s0 t0		# update position of output pointer
	sw t3 0(s4)		# output[t0] = t3
	j Exit

expand_c_bnez:
	li t0 0x380 		# RS1' filter 111|00000|00
	and t2 a2 t0 		# t2 RS1'|00000|00
	li t0 0x400		# ex-RS1 01000|00000|00
	add t2 t2 t0            # RS1|00000|00
	 
	li t0 0x18   		# offset[2|1]filter 11000
	and t1 a2 t0  	        # offset[2|1]0|00

	li t0 0xC00		# offset[4|3]filter 11|000|00000|00
	and t3 a2 t0 		# offset[4|3]|000|00000|00
	slli t1 t1 5            # offset[2|1]0|0000000
	add t1 t1 t3 		# offset[4|1]0|0000000

	li t3 0x63		# opcode 1100011
	add t3 t3 t1 		# offset[4:1]0|opcode

	li t0 0x1000		# func3 001|00000|0000000
	add t3 t3 t0          	# func3|offset[4:1]0|opcode

	slli t2 t2 8            # rs1|000|00000|0000000
	add t3 t3 t2 		# x0|rs1|func3|func3|offset[4:1]0|opcode

	li t0 0x4  		# offset[5]filter 100
	and t1 a2 t0 		# offset[5]|00
	slli t1 t1 23           # offset[5]|00000|00000|000|00000|0000000
	add t3 t3 t1 		# offset[5]|x0|rs1|func3|func3|offset[4:1]0|opcode

	li t0 0x60 		# offset[7|6]filter 1100000
	and t1 a2 t0            # offset[7|6]000|00
	slli t1 t1 21           # offset[7|6]0|00000|00000|000|00000|0000000
	add t3 t3 t1 	 	# offset[7:5]|x0|rs1|func3|func3|offset[4:1]0|opcode

	li t0 0x1000 		# offset[8]filter 100|000|00000|00
	and t1 a2 t0            # offset[8]00|000|00000|00
	srli t1 t1 5
	add t3 t3 t1
	slli t1 t1 24           # offset[8]|0X31
	srai t1 t1 3
	add t3 t3 t1 	 	# offset[8]...offset[7:5]|x0|rs1|func3|func3|offset[4:1]0|opcode


	slli t0 a3 2 		# update the offset
	li s4 0			# s4 is the iterator of output
	add s4 s0 t0		# update position of output pointer
	sw t3 0(s4)		# output[t0] = t3
	j Exit


expand_c_srli:
	li t0 0x7c		# shamt[4:0] filter 11111|00
	and t1 a2 t0 		# shamt[4:0]|00

	li t0 0x380		# RD' filter 111|00000|00
	and t2 a2 t0 		# rd'|00000|00
	li t0 0x400		# ex-RD 01000|00000|00
	add t2 t2 t0            # RD|00000|00

	li t3 0x13		# opcode 0010011
	add t3 t3 t2 		# rd|opcode

	li t0 0x5000		# 101|00000|0000000
	add t3 t3 t0		# func3|rd|opcode

	slli t2 t2 8		# rd|000|00000|0000000
	add t3 t3 t2 		# rd|func3|rd|opcode

	slli t1 t1 18 		# shamt[4:0]|00000|000|00000|0000000
	add t3 t3 t1		# func7|shamt[4:0]|rd|func3|rd|opcode

	slli t0 a3 2 		# update the offset
	li s4 0			# s4 is the iterator of output
	add s4 s0 t0		# update position of output pointer
	sw t3 0(s4)		# output[t0] = t3
	j Exit

expand_c_srai:
	li t0 0x7c		# shamt[4:0] filter 11111|00
	and t1 a2 t0 		# shamt[4:0]|00

	li t0 0x380		# RD' filter 111|00000|00
	and t2 a2 t0 		# rd'|00000|00
	li t0 0x400		# ex-RD 01000|00000|00
	add t2 t2 t0            # RD|00000|00

	li t3 0x13		# opcode 0010011
	add t3 t3 t2 		# rd|opcode

	li t0 0x5000		# 101|00000|0000000
	add t3 t3 t0		# func3|rd|opcode

	slli t2 t2 8		# rd|000|00000|0000000
	add t3 t3 t2 		# rd|func3|rd|opcode

	slli t1 t1 18 		# shamt[4:0]|00000|000|00000|0000000
	add t3 t3 t1		# shamt[4:0]|rd|func3|rd|opcode

	li t0 0x40000000	# func7 0100000|00000|00000|000|00000|0000000
	add t3 t3 t0 		# func7|shamt[4:0]|rd|func3|rd|opcode

	slli t0 a3 2 		# update the offset
	li s4 0			# s4 is the iterator of output
	add s4 s0 t0		# update position of output pointer
	sw t3 0(s4)		# output[t0] = t3
	j Exit


expand_c_andi:
	li t0 0x7c		# imm[4:0] filter 11111|00
	and t1 a2 t0 		# imm[4:0]|00

	li t0 0x380		# RD' filter 111|00000|00
	and t2 a2 t0 		# rd'|00000|00
	li t0 0x400		# ex-RD 01000|00000|00
	add t2 t2 t0            # RD|00000|00

	li t3 0x13		# opcode 0010011
	add t3 t3 t2 		# rd|opcode

	li t0 0x7000		# 111|00000|0000000
	add t3 t3 t0		# func3|rd|opcode

	slli t2 t2 8		# rd|000|00000|0000000
	add t3 t3 t2 		# rd|func3|rd|opcode

	slli t1 t1 18 		# imm[4:0]|00000|000|00000|0000000
	add t3 t3 t1		# imm[4:0]|rd|func3|rd|opcode

	li t0 0x1000		# imm[5] filter 1|00|000|00000|00
	and t0 a2 t0
	slli t0 t0 19           # imm[5]|31X0
	srai t0 t0 6
	add t3 t3 t0 		# imm[5]...imm[4:0]|rd|func3|rd|opcode

	slli t0 a3 2 		# update the offset
	li s4 0			# s4 is the iterator of output
	add s4 s0 t0		# update position of output pointer
	sw t3 0(s4)		# output[t0] = t3
	j Exit


expand_c_j:
	li t0 0x6f		# opcode 1101111
	
	j expand_c_j_jal_imm

expand_c_jal:
	li t0 0x6f		# opcode 1101111
	li t1 0x80		# rd is 00001|0000000
	add t0 t0 t1		# t0 is rd|opcode
	
	j expand_c_j_jal_imm

expand_c_j_jal_imm:
	li t1 0x38		# mask of imm[3:1] --> 0011 1000
	and t2 a2 t1		# t2 is imm[3:1] >> 18
	slli t2 t2 18		# t2 is imm[3:1]
	add t0 t0 t2		# finish imm[3:1]

	li t1 0x800		# mask of imm[4] --> 1000 0000 0000
	and t2 a2 t1 		# t2 is imm[4] >> 13
	slli t2 t2 13		# t2 is imm[4]
	add t0 t0 t2		# finish imm[4:1]

	li t1 0x4		# mask of imm[5] --> 0100
	and t2 a2 t1		# t2 is imm[5] >> 23
	slli t2 t2 23		# t2 is imm[5]
	add t0 t0 t2		# finish imm[5]

	li t1 0x80		# mask of imm[6] --> 1000 0000
	and t2 a2 t1 		# t2 is imm[6] >> 19
	slli t2 t2 19		# t2 is imm[6]
	add t0 t0 t2		# finish imm[6]

	li t1 0x40		# mask of imm[7] --> 0100 0000
	and t2 a2 t1 		# t2 is imm[7] >> 21
	slli t2 t2 21		# t2 is imm[7]
	add t0 t0 t2		# finish imm[7]

	li t1 0x600		# mask of imm[9:8] --> 0110 0000 0000
	and t2 a2 t1 		# t2 is imm[9:8] >> 19
	slli t2 t2 19		# t2 is imm[9:8]
	add t0 t0 t2		# finish imm[9:8]

	li t1 0x100		# mask of imm[10] --> 0001 0000 0000
	and t2 a2 t1		# t2 is imm[10] >> 22
	slli t2 t2 22		# t2 is imm[10]
	add t0 t0 t2 		# finish imm[10]

	li t1 0x1000		# mask of imm[11] --> 0001 0000 0000 0000
	and t2 a2 t1		# t2 is imm[11] >> 8
	slli t2 t2 8		# t2 is imm[11]
	add t0 t0 t2		# finish imm[11]

	beq t2 x0 expand_c_j_jal_end
	li t1 0x800ff000	# mask of imm[20|19:12]
	or t0 t0 t1		# sign extention
expand_c_j_jal_end:
	slli t1 a3 2 		# update the offset
	li s4 0			# s4 is the iterator of output
	add s4 s0 t1		# update position of output pointer
	sw t0 0(s4)		# output[t0] = t0
	j Exit


# Define the following coding:
# c.add	--> 1
# c.mv	--> 2
# c.jr	--> 3
# c.jalr--> 4
# c.li	--> 5
# c.lui	--> 6
# c.addi--> 7
# c.slli--> 8
# c.lw	--> 9
# c.sw	--> 10
# c.and	--> 11
# c.or	--> 12
# c.xor	--> 13
# c.sub	--> 14
# c.beqz--> 15
# c.bneq--> 16
# c.srli--> 17
# c.srai--> 18
# c.andi--> 19
# c.j	--> 20
# c.jal	--> 21