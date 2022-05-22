#include <stdio.h>
#include <string.h>

#include "utils.h"

/* Your code here... */

/*read data line by line*/
int read_line(FILE *in, int *out) {
    int num = fgetc(in); /*get number*/
    int i = 31; /*iterator*/
    int j; /*iterator*/
    while(1) {
        if (num == '1') out[i--] = 1; /*when num is 1*/
        if (num == '0') out[i--] = 0; /*when num is 0*/
        if (num == EOF) { /*when touch the end of the file*/
            if (out[0] == -1) {
                for (j = 0; j < 32; ++j) {
                    out[j] = -1; /*initialize out*/
                }
            }
            return 0;
        }
        if (num == '\n') return 1; /*when touch the end of a line*/
        num = fgetc(in);
    }
    /*int i;
    int tmp;
    *//*edge case*//*
    if (!in || !out) {
        printf("the input file or the receiver array is empty\n");
        return 0;
    }
    *//*reverse the input data:
     * let the last significant bit equal to out[0]
     * and the most significant bit be out[31]*//*

    for (i = 31; i >= 0; i--) {
        tmp = fgetc(in); *//*temp result*//*
        if (tmp == EOF) {*//*if touch the end of the file*//*
            memset(out, 0, 32 * sizeof(int));
            return 0;
        }
        out[i] = tmp; *//*read in data*//*
        out[i] -= 48;
    }
    *//*use the additional "fgetc" to get the '\n' for each line*//*
    tmp = fgetc(in);
    while (tmp == '\t') { *//* deal with '\t'*//*
        tmp = fgetc(in);
        if (tmp == '\n') break;
    }
    return 1;*/
}

/*write data line by line*/
int write_line(FILE *out, const int *in, int mode) {
    /*Determine the value of i via the value fo mode.*/
    int i = (mode == 0) ? 31 : 15;
    /*Deal with edge cases*/
    if (!out || !in) {
        printf("the output file or the input array is empty\n");
        return 1;
    }
    /*Write data into file "out" reversely*/
    for (; i > 0; --i) {
        fprintf(out, "%d", in[i]);
    }
    /*the last element of the line should be followed by '\n'*/
    fprintf(out, "%d\n", in[0]);
    return 0;
}

/*turn binary number into decimal number*/
int get_Decimal_number(const int *data, int low, int high) {
    int i, res = 0;/*iterator and result*/
    for (i = low; i <= high; ++i) {
        res += (data[i] << (i - low));
    }
    return res;
}

/*get s type imm*/
int get_S_type_imm(instruction *line) {
    if (line->data[31] == 0) {
        return line->rd + (get_Decimal_number(line->data, 25, 30) << 5);
        /*get imm[10:5] + imm[4:0]*/
    } else {
        long imm = line->rd + (line->funct7 << 5);
        /*get imm[11:5] + imm[4:0]*/
        imm |= 0xfffff000;
        /* 1111 1111 1111 1111 1111 0000 0000 0000*/
        return (int) imm;
    }
}

/*get sb type imm*/
int get_SB_type_imm(instruction *line) {
    long imm = get_Decimal_number(line->data, 8, 11);
    /*get imm[4:1]*/
    imm += (get_Decimal_number(line->data, 25, 30) << 4);
    /*get imm[10:5] + imm[4:1]*/
    imm += (line->data[7] << 10);
    /*get imm[11] + imm[10:5] + imm[4:1]*/
    if (line->data[31] == 0) {
        return (int) imm;
    } else {
        imm += (line->data[31] << 11);
        /*get imm[12] + imm[11] + imm[10:5] + imm[4:1]*/
        imm |= 0xfffff000;
        /*1111 1111 1111 1111 1111 0000 0000 0000*/
        return (int) imm;
    }
}

/*get uj type imm*/
int get_UJ_type_imm(instruction *line) {
    long imm = get_Decimal_number(line->data, 21, 30);
    /*get imm[10:1]*/
    imm += (line->data[20] << 10);
    /*get imm[11] + imm[10:1]*/
    imm += (get_Decimal_number(line->data, 12, 19) << 11);
    /*get imm[19:12] + imm[11] + imm[10:1]*/
    if (line->data[31] == 0) {
        return (int) imm;
    } else {
        imm += (line->data[31] << 19);
        /*get imm[20] + imm[19:12] + imm[11] + imm[10:1]*/
        imm |= 0xfff00000;
        /* 1111 1111 1111 0000 0000 0000 0000 0000*/
        return (int) imm;
    }
}

/*get i type unsigned imm*/
/*int get_I_type_unsigned_imm(instruction *line) {
    return get_Decimal_number(line->data, 20, 31);
}*/

/*get i type imm*/
int get_I_type_imm(instruction *line) {
    if (line->data[31] == 0) {
        return get_Decimal_number(line->data, 20, 30);
        /*get imm[10:0]*/
    } else {
        long imm = get_Decimal_number(line->data, 20, 31);
        /*get imm[11:0]*/
        imm |= 0xfffff000;
        /*1111 1111 1111 1111 1111 0000 0000 0000*/
        return (int) imm;
    }
}

/*get u type imm*/
int get_U_type_imm(instruction *line) {
    if (line->data[31] == 0) {
        return get_Decimal_number(line->data, 12, 30);
        /*get imm[10:0]*/
    } else {
        long imm = get_Decimal_number(line->data, 12, 31);
        /*get imm[11:0]*/
        imm |= 0xfff00000;
        /*1111 1111 1111 0000 0000 0000 0000 0000*/
        return (int) imm;
    }
}

/*check if reg is one of x8/x9/x10/x11/x12/x13/x14/x15*/
int check_reg(int reg) {
    return (reg == 8 || reg == 9 || reg == 10 || reg == 11 || reg == 12 || reg == 13 || reg == 14 || reg == 15);
}


