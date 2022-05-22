.section .text
.globl draw
.type draw, @function
draw:
	li a0, 0
	li a1, 40
	li a2, 160
	li a3, 40
	li a4, 0xEA69
	jal LCD_DrawLine
	li a0, 80
	li a1, 20
	li a2, 20
	li a3, 0x074E
	jal LCD_DrawCircle
	li a0, 80
	li a1, 60
	li a2, 2
	li a3, 16
	li a4, 0xEDC0
	jal LCD_ShowChinese
	li a0, 10
	ecall
	# YOUR CODE HERE
