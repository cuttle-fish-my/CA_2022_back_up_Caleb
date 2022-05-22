#ifndef UTILS_H
#define UTILS_H

#include "stdio.h"
#include "compression.h"

/* Your code here... */




int read_line(FILE *in, int *out);

/*read_line is used for reading data from file "in" line by line
 * the int pointer "out" is a pointer to an int array
 * with length 32 to store data from a line*/

int write_line(FILE *out, const int *in, int mode);

/*write_line is used for writing data from "in" to "out" line by line
 * "in" is a pointer to an int array with length 32.
 * The parameter "mode" is used to determined whether the instruction is compressed or not
 * mode == 0 means this line is not compressed(32-bit long) and mode == 1 means the line is 16-bit long*/

int get_Decimal_number(const int *data, int low, int high);

/*turn binary number into decimal number*/

int get_S_type_imm(instruction *line);

/*get s type imm*/

int get_SB_type_imm(instruction *line);

/*get sb type imm*/

int get_UJ_type_imm(instruction *line);

/*get uj type imm*/

int get_I_type_imm(instruction *line);

/*get i type imm*/

int get_I_type_unsigned_imm(instruction *line);

/*get i type unsigned imm*/

int get_U_type_imm(instruction *line);

/*get u type unsigned imm*/

int check_reg(int reg);

/*this func is used to check whether the register is one of x8 to x15*/
#endif