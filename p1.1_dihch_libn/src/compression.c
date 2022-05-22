#include <string.h>

#include "compression.h"
#include "utils.h"

/* Your code here... */
void check_type(instruction *line) {
    line->opcode = get_Decimal_number(line->data, 0, 6); /*opcode*/
    line->funct3 = get_Decimal_number(line->data, 12, 14); /*funct3*/
    line->funct7 = get_Decimal_number(line->data, 25, 31); /*funct7*/
    line->rd = get_Decimal_number(line->data, 7, 11); /*rd*/
    line->rs1 = get_Decimal_number(line->data, 15, 19); /*rs1*/
    line->rs2 = get_Decimal_number(line->data, 20, 24); /*rs2*/
    line->shamt = get_Decimal_number(line->data, 20, 24); /*shamt*/
    switch (line->opcode) {
        case 0x33: /*R-type*/
            /*add / and / or / xor / sub*/
            switch (line->funct3) {
                case 0x0:
                    /*add / sub*/
                    switch (line->funct7) {
                        case 0x0:
                            /*add*/
                            check_add(line, line->rd, line->rs1, line->rs2);
                            break;
                        case 0x20:
                            /*sub*/
                            check_sub(line, line->rd, line->rs1, line->rs2);
                            break;
                        default:/*other that can not be compressed*/
                            line->type = none;
                            break;
                    }
                    break;
                case 0x4:
                    /*xor*/
                    if (line->funct7 == 0x0) check_xor(line, line->rd, line->rs1, line->rs2);
                    else line->type = none;
                    break;
                case 0x6:
                    /*or*/
                    if (line->funct7 == 0x0) check_or(line, line->rd, line->rs1, line->rs2);
                    else line->type = none;
                    break;
                case 0x7:
                    /*and*/
                    if (line->funct7 == 0x0) check_and(line, line->rd, line->rs1, line->rs2);
                    else line->type = none;
                    break;
                default:
                    /*other that can not be compressed*/
                    line->type = none;
                    break;
            }
            break;
        case 0x37: /*U-type*/
            /*lui*/
            line->imm = get_U_type_imm(line);
            check_lui(line, line->rd, line->imm); /*imm is signed*/
            break;
        case 0x13: /*I-type*/
            /*slli / addi / srli / srai / andi*/
            line->imm = get_I_type_imm(line); /*get imm for all I type instructions below */
            switch (line->funct3) {
                case 0x0:
                    /*addi*/
                    check_addi(line, line->rd, line->rs1, line->imm); /*imm is signed*/
                    break;
                case 0x1:
                    /*slli*/
                    if (line->funct7 == 0x0) check_slli(line, line->rd, line->rs1); /*shamt is unsigned*/
                    else line->type = none;
                    break;
                case 0x5:
                    /*srli / srai*/
                    switch (line->funct7) {
                        case 0x0:
                            /*srli*/
                            check_srli(line, line->rd, line->rs1);
                            break;
                        case 0x20:
                            /*srai*/
                            check_srai(line, line->rd, line->rs1);
                            break;
                        default: /*other that can not be compressed*/
                            line->type = none;
                            break;
                    }
                    break;
                case 0x7:
                    /*andi*/
                    check_andi(line, line->rd, line->rs1, line->imm);
                    break;
                default:
                    /*other that can not be compressed*/
                    line->type = none;
                    break;
            }
            break;
        case 0x6f: /*UJ-type*/
            /*jal*/
            line->imm = get_UJ_type_imm(line);
            check_jal(line, line->rd);
            /*no need to check whether imm fits in 6-bits binary number*/
            break;
        case 0x67: /*I-type*/
            /*jalr*/
            line->imm = get_I_type_imm(line);
            check_jalr(line, line->rd, line->rs1, line->imm);
            break;
        case 0x03: /*I-type*/
            /*lw*/
            line->imm = get_I_type_imm(line);
            if (line->funct3 == 0x2) { /*funct3*/
                check_lw(line, line->rd, line->rs1, line->imm);
            } else line->type = none;
            break;
        case 0x23: /*S-type*/
            /*sw*/
            line->imm = get_S_type_imm(line);
            if (line->funct3 == 0x2) { /*funct3*/
                check_sw(line, line->rs1, line->imm);
            } else line->type = none;
            break;
        case 0x63: /*SB-type*/
            /*beq / bne*/
            line->imm = get_SB_type_imm(line);
            switch (line->funct3) {
                case 0x0:
                    /*beq*/
                    check_beq(line, line->rs1, line->rs2);
                    break;
                case 0x1:
                    /*bne*/
                    check_bne(line, line->rs1, line->rs2);
                    break;
                case 0x4:
                    /*blt*/
                    line->type = none_blt;
                    break;
                case 0x5:
                    /*bge*/
                    line->type = none_bge;
                    break;
                case 0x6:
                    /*bltu*/
                    line->type = none_bltu;
                    break;
                case 0x7:
                    /*bgeu*/
                    line->type = none_bgeu;
                    break;
                default:
                    /*other that can not be compressed*/
                    line->type = none;
                    break;
            }
            break;
        default:
            /*other that can not be compressed*/
            line->type = none;
            break;
    }
}

/*check add constrains*/
void check_add(instruction *line, int rd, int rs1, int rs2) {
    if ((rd == rs1 && rd != 0 && rs2 != 0) || (rs1 == 0 && rd != 0 && rs2 != 0)) line->type = add;
    else line->type = none;
}

/*check lui constrains*/
void check_lui(instruction *line, int rd, int imm) {
    if (rd != 0 && rd != 2 && imm >= -32 && imm <= 31) line->type = lui;
        /*rd != x0 or x2 and -32 <= imm <=31  (6-bit) */
    else line->type = none;
}

/*check slli constrains*/
void check_slli(instruction *line, int rd, int rs1) {
    if (rd == rs1 && rd != 0) line->type = slli;
    else line->type = none;
}

/*check addi constrains*/
void check_addi(instruction *line, int rd, int rs1, int imm) {
    if ((rs1 == 0 && rd != 0 && imm >= -32 && imm <= 31) ||
        (rd == rs1 && rd != 0 && imm != 0 && imm >= -32 && imm <= 31))
        line->type = addi;
        /*if rs1 = x0 and rd != x0 and imm fits in 6-bit signed number, then it is c.li format*/
        /*if rd = rs1 and rd != x0 and imm fits in 6-bit signed number and imm != 0, then it is c.addi*/
    else line->type = none;
}

/*check srli constrains*/
void check_srli(instruction *line, int rd, int rs1) {
    if (rd == rs1 && check_reg(rd) && check_reg(rs1)) line->type = srli;
    else line->type = none;
}

/*check srai constrains*/
void check_srai(instruction *line, int rd, int rs1) {
    if (rd == rs1 && check_reg(rd) && check_reg(rs1)) line->type = srai;
    else line->type = none;
}

/*check jalr constrains*/
void check_jalr(instruction *line, int rd, int rs1, int imm) {
    if ((rd == 0 || rd == 1) && rs1 != 0 && imm == 0) line->type = jalr;
        /*if rd = x0 and rs1 != x0 and imm = 0, then it is c.jr format*/
        /*if rd = x1 and rs1 != x0 and imm = 0, then it is c.jalr format*/
    else line->type = none;
}

/*check andi constrains*/
void check_andi(instruction *line, int rd, int rs1, int imm) {
    if (imm >= -32 && imm <= 31 && rd == rs1 && check_reg(rd) && check_reg(rs1)) line->type = andi;
        /*if imm fits in 6-bits binary number, then it is c.andi format*/
    else line->type = none;
}

/*check lw constrains*/
void check_lw(instruction *line, int rd, int rs1, int imm) {
    if (check_reg(rd) && check_reg(rs1) && imm >= 0 && imm <= 124 && (imm % 4) == 0) line->type = lw;
    /*if imm is positive and imm is a multiple of 4 and 0 <= imm <= 127, then it is lw type*/
    else line->type = none;
    /*else it is incompressible*/
}

/*check sw constrains*/
void check_sw(instruction *line, int rs1, int imm) {
    if (check_reg(rs1) && check_reg(line->rs2) && imm >= 0 && imm <= 124 && (imm % 4) == 0) line->type = sw;
    /*if imm is positive and imm is a multiple of 4 and 0 <= imm <= 127, then it is lw type*/
    else line->type = none;
    /*else it is incompressible*/
}

/*check beq constrains*/
void check_beq(instruction *line, int rs1, int rs2) {
    if (check_reg(rs1) && rs2 == 0) line->type = beq;
    else line->type = none_beq;
}

/*check bne constrains*/
void check_bne(instruction *line, int rs1, int rs2) {
    if (check_reg(rs1) && rs2 == 0) line->type = bne;
    else line->type = none_bne;
}

/*check and constrains*/
void check_and(instruction *line, int rd, int rs1, int rs2) {
    if (rs1 == rd && check_reg(rd) && check_reg(rs1) && check_reg(rs2)) line->type = and;
    else line->type = none;
}

/*check or constrains*/
void check_or(instruction *line, int rd, int rs1, int rs2) {
    if (rs1 == rd && check_reg(rd) && check_reg(rs1) && check_reg(rs2)) line->type = or;
    else line->type = none;
}

/*check xor constrains*/
void check_xor(instruction *line, int rd, int rs1, int rs2) {
    if (rs1 == rd && check_reg(rd) && check_reg(rs1) && check_reg(rs2)) line->type = xor;
    else line->type = none;
}

/*check sub constrains*/
void check_sub(instruction *line, int rd, int rs1, int rs2) {
    if (rs1 == rd && check_reg(rd) && check_reg(rs1) && check_reg(rs2)) line->type = sub;
    else line->type = none;
}

/*check jal constrains*/
void check_jal(instruction *line, int rd) {
    if (rd == 0 || rd == 1) line->type = jal;
    else line->type = none_jal;
}

/*compress CR format*/
compressed_instruction *compression_CR(instruction *line, compressed_instruction *_line) {
    _line->type = line->type;
    _line->data[0] = 0;/*opcode*/
    _line->data[1] = 1;/*opcode*/
    if (_line->type == add)/*C.ADD or C.MV*/
    {
        _line->data[15] = 1;/*FUNCT4*/
        _line->data[14] = 0;/*FUNCT4*/
        _line->data[13] = 0;/*FUNCT4*/
        _line->data[12] = 1;/*FUNCT4*/
        if (line->rs1 == 0)/*C.MOV*/
        {
            _line->data[12] = 0;/*FUNCT4*/
        }

        _line->data[11] = line->data[11];/*RD*/
        _line->data[10] = line->data[10];/*RD*/
        _line->data[9] = line->data[9];/*RD*/
        _line->data[8] = line->data[8];/*RD*/
        _line->data[7] = line->data[7];/*RD*/

        _line->data[6] = line->data[24];/*RS2*/
        _line->data[5] = line->data[23];/*RS2*/
        _line->data[4] = line->data[22];/*RS2*/
        _line->data[3] = line->data[21];/*RS2*/
        _line->data[2] = line->data[20];/*RS2*/
        return _line;
    }
    if (_line->type == jalr)/*C.JR or C.JALR*/
    {
        _line->data[15] = 1;/*FUNCT4*/
        _line->data[14] = 0;/*FUNCT4*/
        _line->data[13] = 0;/*FUNCT4*/
        _line->data[12] = 1;/*FUNCT4*/
        if (line->rd == 0)/*C.JR*/
        {
            _line->data[12] = 0;/*FUNCT4*/
        }

        _line->data[11] = line->data[19];/*RS1*/
        _line->data[10] = line->data[18];/*RS1*/
        _line->data[9] = line->data[17];/*RS1*/
        _line->data[8] = line->data[16];/*RS1*/
        _line->data[7] = line->data[15];/*RS1*/

        _line->data[6] = 0;/*RS2*/
        _line->data[5] = 0;/*RS2*/
        _line->data[4] = 0;/*RS2*/
        _line->data[3] = 0;/*RS2*/
        _line->data[2] = 0;/*RS2*/
        return _line;
    }
    return NULL;
}

/*compress CI format*/
compressed_instruction *compression_CI(instruction *line, compressed_instruction *_line) {
    _line->type = line->type;
    _line->data[11] = line->data[11];/*RD*/
    _line->data[10] = line->data[10];/*RD*/
    _line->data[9] = line->data[9];/*RD*/
    _line->data[8] = line->data[8];/*RD*/
    _line->data[7] = line->data[7];/*RD*/


    if (_line->type == addi && line->rs1 == 0)/*C.LI*/
    {
        _line->data[15] = 0;/*FUNCT3*/
        _line->data[14] = 1;/*FUNCT3*/
        _line->data[13] = 0;/*FUNCT3*/
        _line->data[12] = line->data[31];/*imm[5]*/
        _line->data[6] = line->data[24];/*imm[4:0]*/
        _line->data[5] = line->data[23];/*imm[4:0]*/
        _line->data[4] = line->data[22];/*imm[4:0]*/
        _line->data[3] = line->data[21];/*imm[4:0]*/
        _line->data[2] = line->data[20];/*imm[4:0]*/
        _line->data[0] = 1;/*opcode*/
        _line->data[1] = 0;/*opcode*/
        return _line;
    }

    if (_line->type == lui)/*C.LUI*/
    {
        _line->data[15] = 0;/*FUNCT3*/
        _line->data[14] = 1;/*FUNCT3*/
        _line->data[13] = 1;/*FUNCT3*/
        _line->data[12] = line->data[31];/*nzimm[17]*/
        _line->data[6] = line->data[16];/*nzimm[16:12]*/
        _line->data[5] = line->data[15];/*nzimm[16:12]*/
        _line->data[4] = line->data[14];/*nzimm[16:12]*/
        _line->data[3] = line->data[13];/*nzimm[16:12]*/
        _line->data[2] = line->data[12];/*nzimm[16:12]*/
        _line->data[0] = 1;/*opcode*/
        _line->data[1] = 0;/*opcode*/
        return _line;
    }

    if (_line->type == addi && line->rd == line->rs1)/*C.ADDI*/
    {
        _line->data[15] = 0;/*FUNCT3*/
        _line->data[14] = 0;/*FUNCT3*/
        _line->data[13] = 0;/*FUNCT3*/
        _line->data[12] = line->data[31];/*imm[5]*/
        _line->data[6] = line->data[24];/*imm[4:0]*/
        _line->data[5] = line->data[23];/*imm[4:0]*/
        _line->data[4] = line->data[22];/*imm[4:0]*/
        _line->data[3] = line->data[21];/*imm[4:0]*/
        _line->data[2] = line->data[20];/*imm[4:0]*/
        _line->data[0] = 1;/*opcode*/
        _line->data[1] = 0;/*opcode*/
        return _line;
    }

    if (_line->type == slli)/*C.SLLI*/
    {
        _line->data[15] = 0;/*FUNCT3*/
        _line->data[14] = 0;/*FUNCT3*/
        _line->data[13] = 0;/*FUNCT3*/
        _line->data[12] = 0;/*shamt[5]*/ /*just set 0*/
        _line->data[6] = line->data[24];/*shamt[4:0]*/
        _line->data[5] = line->data[23];/*shamt[4:0]*/
        _line->data[4] = line->data[22];/*shamt[4:0]*/
        _line->data[3] = line->data[21];/*shamt[4:0]*/
        _line->data[2] = line->data[20];/*shamt[4:0]*/
        _line->data[0] = 0;/*opcode*/
        _line->data[1] = 1;/*opcode*/
        return _line;
    }
    return NULL;
}

/*compress CL format*/
compressed_instruction *compression_CL(instruction *line, compressed_instruction *_line) {
    _line->type = line->type;
    _line->data[15] = 0;/*FUNCT3*/
    _line->data[14] = 1;/*FUNCT3*/
    _line->data[13] = 0;/*FUNCT3*/
    _line->data[12] = line->data[25];/*offset[5:3]*/
    _line->data[11] = line->data[24];/*offset[5:3]*/
    _line->data[10] = line->data[23];/*offset[5:3]*/
    _line->data[9] = line->data[17];/*RS1'*/
    _line->data[8] = line->data[16];/*RS1'*/
    _line->data[7] = line->data[15];/*RS1'*/
    _line->data[6] = line->data[22];/*offset[2|6]*/
    _line->data[5] = line->data[26];/*offset[2|6]*/
    _line->data[4] = line->data[9];/*RD'*/
    _line->data[3] = line->data[8];/*RD'*/
    _line->data[2] = line->data[7];/*RD'*/
    _line->data[0] = 0;/*opcode*/
    _line->data[1] = 0;/*opcode*/
    return _line;
}

/*compress CS format*/
compressed_instruction *compression_CS(instruction *line, compressed_instruction *_line) {
    _line->type = line->type;
    /*CS-TYPE1*/
    if (_line->type == sw) {
        _line->data[15] = 1;/*FUNCT3*/
        _line->data[14] = 1;/*FUNCT3*/
        _line->data[13] = 0;/*FUNCT3*/
        _line->data[12] = line->data[25];/*offset[5:3]*/
        _line->data[11] = line->data[11];/*offset[5:3]*/
        _line->data[10] = line->data[10];/*offset[5:3]*/
        _line->data[9] = line->data[17];/*RS1'*/
        _line->data[8] = line->data[16];/*RS1'*/
        _line->data[7] = line->data[15];/*RS1'*/
        _line->data[6] = line->data[9];/*offset[2|6]*/
        _line->data[5] = line->data[26];/*offset[2|6]*/
        _line->data[4] = line->data[22];/*RS2'*/
        _line->data[3] = line->data[21];/*RS2'*/
        _line->data[2] = line->data[20];/*RS2'*/
        _line->data[0] = 0;/*opcode*/
        _line->data[1] = 0;/*opcode*/
        return _line;
    }
    /*CS-TYPE2*/
    _line->data[15] = 1;/*FUNCT6*/
    _line->data[14] = 0;/*FUNCT6*/
    _line->data[13] = 0;/*FUNCT6*/
    _line->data[12] = 0;/*FUNCT6*/
    _line->data[11] = 1;/*FUNCT6*/
    _line->data[10] = 1;/*FUNCT6*/
    _line->data[9] = line->data[17];/*RD'*/
    _line->data[8] = line->data[16];/*RD'*/
    _line->data[7] = line->data[15];/*RD'*/
    if (_line->type == and) {
        _line->data[6] = 1;/*FNCT2*/
        _line->data[5] = 1;/*FNCT2*/
    }
    if (_line->type == or) {
        _line->data[6] = 1;/*FNCT2*/
        _line->data[5] = 0;/*FNCT2*/
    }
    if (_line->type == xor) {
        _line->data[6] = 0;/*FNCT2*/
        _line->data[5] = 1;/*FNCT2*/
    }
    if (_line->type == sub) {
        _line->data[6] = 0;/*FNCT2*/
        _line->data[5] = 0;/*FNCT2*/
    }
    _line->data[4] = line->data[22];/*RS2'*/
    _line->data[3] = line->data[21];/*RS2'*/
    _line->data[2] = line->data[20];/*RS2'*/
    _line->data[0] = 1;/*opcode*/
    _line->data[1] = 0;/*opcode*/
    return _line;
}

/*compress CB format*/
compressed_instruction *
compression_CB(instruction *line, const int *address, int _line_No, compressed_instruction *_line) {
    _line->type = line->type;
    _line->data[0] = 1;/*opcode*/
    _line->data[1] = 0;/*opcode*/
    /*CB-TYPE1*/
    if (_line->type == beq || _line->type == bne) {
        int imm = line->imm;
        int new_imm = address[_line_No + imm / 2] - address[_line_No];
        int i = 1, shift = 0;
        _line->data[15] = 1;/*FUNCT3*/
        _line->data[14] = 1;/*FUNCT3*/
        _line->data[13] = 0;/*FUNCT3*/
        if (_line->type == bne) _line->data[13] = 1;

        /*IMM*/
        _line->data[3] = ((new_imm & i) >> (shift++)) & 1;
        i = i << 1;/*#1*/
        _line->data[4] = ((new_imm & i) >> (shift++)) & 1;
        i = i << 1;/*#2*/
        _line->data[10] = ((new_imm & i) >> (shift++)) & 1;
        i = i << 1;/*#3*/
        _line->data[11] = ((new_imm & i) >> (shift++)) & 1;
        i = i << 1;/*#4*/
        _line->data[2] = ((new_imm & i) >> (shift++)) & 1;
        i = i << 1;/*#5*/
        _line->data[5] = ((new_imm & i) >> (shift++)) & 1;
        i = i << 1;/*#6*/
        _line->data[6] = ((new_imm & i) >> shift) & 1;/*#7*/
        _line->data[12] = ((new_imm & (1 << 31)) >> 31) & 1;/*#8*/

        /*RD'/RS1'*/
        _line->data[7] = line->data[15];/*RD'/RS1'*/
        _line->data[8] = line->data[16];/*RD'/RS1'*/
        _line->data[9] = line->data[17];/*RD'/RS1'*/
        return _line;
    }
    /*CB-TYPE2*/
    if (_line->type == srli || _line->type == srai || _line->type == andi) {
        _line->data[15] = 1;/*FUNCT3*/
        _line->data[14] = 0;/*FUNCT3*/
        _line->data[13] = 0;/*FUNCT3*/

        /*FUNCT2*/
        if (_line->type == srli) {
            _line->data[11] = 0;/*FUNCT2*/
            _line->data[10] = 0;/*FUNCT2*/
        }
        /*FUNCT2*/
        if (_line->type == srai) {
            _line->data[11] = 0;/*FUNCT2*/
            _line->data[10] = 1;/*FUNCT2*/
        }
        /*FUNCT2*/
        if (_line->type == andi) {
            _line->data[11] = 1;/*FUNCT2*/
            _line->data[10] = 0;/*FUNCT2*/
        }

        /*IMM*/ /*shamt*/
        _line->data[2] = line->data[20];/*#0*/
        _line->data[3] = line->data[21];/*#1*/
        _line->data[4] = line->data[22];/*#2*/
        _line->data[5] = line->data[23];/*#3*/
        _line->data[6] = line->data[24];/*#4*/
        _line->data[12] = 0;
        if (_line->type == andi) {
            _line->data[12] = line->data[31];/*#5*/
        }

        /*RD'/RS1'*/
        _line->data[7] = line->data[7];/*RD'/RS1'*/
        _line->data[8] = line->data[8];/*RD'/RS1'*/
        _line->data[9] = line->data[9];/*RD'/RS1'*/
        return _line;
    }
    return NULL;
}

/*compress CJ*/
compressed_instruction *
compression_CJ(instruction *line, const int *address, int _line_No, compressed_instruction *_line) {
    int imm = line->imm; /*imm*/
    int new_imm = address[_line_No + imm / 2] - address[_line_No]; /*new_imm*/
    int i = 1, shift = 0;
    _line->type = line->type;
    _line->data[0] = 1;/*OPCODE*/
    _line->data[1] = 0;/*OPCODE*/
    /*FUNCT3*/
    if (line->rd == 0) /*C.J*/
    {
        _line->data[15] = 1;/*FUNCT3*/
        _line->data[14] = 0;/*FUNCT3*/
        _line->data[13] = 1;/*FUNCT3*/
    } else /*C.JAL*/
    {
        _line->data[15] = 0;/*FUNCT3*/
        _line->data[14] = 0;/*FUNCT3*/
        _line->data[13] = 1;/*FUNCT3*/
    }

    /*JUMP TARGET*/
    _line->data[3] = ((new_imm & i) >> (shift++)) & 1;
    i = i << 1;/*#1*/
    _line->data[4] = ((new_imm & i) >> (shift++)) & 1;
    i = i << 1; /*#2*/
    _line->data[5] = ((new_imm & i) >> (shift++)) & 1;
    i = i << 1;/*#3*/
    _line->data[11] = ((new_imm & i) >> (shift++)) & 1;
    i = i << 1;/*#4*/
    _line->data[2] = ((new_imm & i) >> (shift++)) & 1;
    i = i << 1;/*#5*/
    _line->data[7] = ((new_imm & i) >> (shift++)) & 1;
    i = i << 1;/*#6*/
    _line->data[6] = ((new_imm & i) >> (shift++)) & 1;
    i = i << 1;/*#7*/
    _line->data[9] = ((new_imm & i) >> (shift++)) & 1;
    i = i << 1;/*#8*/
    _line->data[10] = ((new_imm & i) >> (shift++)) & 1;
    i = i << 1;/*#9*/
    _line->data[8] = ((new_imm & i) >> shift) & 1;/*#10*/
    _line->data[12] = ((new_imm & (1 << 31)) >> 31) & 1;/*#11*/

    return _line;
}

/*change offsets of bltu/bgeu/bge/blt*/
void change_offset(instruction *line, const int *address, int _line_No) {
    int imm = line->imm; /*imm*/
    int new_imm = address[_line_No + imm / 2] - address[_line_No]; /*new_imm*/
    int i = 1, shift = 0;
    line->data[8] = ((new_imm & i) >> (shift++)) & 1;/*#1*/
    i = i << 1;/*#1*/
    line->data[9] = ((new_imm & i) >> (shift++)) & 1;/*#2*/
    i = i << 1;/*#2*/
    line->data[10] = ((new_imm & i) >> (shift++)) & 1;/*#3*/
    i = i << 1;/*#3*/
    line->data[11] = ((new_imm & i) >> (shift++)) & 1;/*#4*/
    i = i << 1;/*#4*/
    line->data[25] = ((new_imm & i) >> (shift++)) & 1;/*#5*/
    i = i << 1;/*#5*/
    line->data[26] = ((new_imm & i) >> (shift++)) & 1;/*#6*/
    i = i << 1;/*#6*/
    line->data[27] = ((new_imm & i) >> (shift++)) & 1;/*#7*/
    i = i << 1;/*#7*/
    line->data[28] = ((new_imm & i) >> (shift++)) & 1;/*#8*/
    i = i << 1;/*#8*/
    line->data[29] = ((new_imm & i) >> (shift++)) & 1;/*#9*/
    i = i << 1;/*#9*/
    line->data[30] = ((new_imm & i) >> (shift++)) & 1;/*#10*/
    i = i << 1;/*#10*/
    line->data[7] = ((new_imm & i) >> (shift++)) & 1;/*#11*/
    i = i << 1;/*#11*/

    line->data[31] = ((new_imm & (1 << 31)) >> 31) & 1;/*#12*/

/*    if (line->type == none_bltu || line->type == none_bgeu) {
        line->data[31] = ((new_imm & i) >> shift) & 1;*//*#12*//*
    }*/
}

/*void change_offset_jal(instruction *line, const int *address, int _line_No) {
    int imm = line->imm;
    int new_imm = address[_line_No + imm / 2] - address[_line_No];
    int i = 1, shift = 0;
    line->data[21] = ((new_imm & i) >> (shift++)) & 1;*//*#1*//*
    i = i << 1;
    line->data[22] = ((new_imm & i) >> (shift++)) & 1;*//*#2*//*
    i = i << 1;
    line->data[23] = ((new_imm & i) >> (shift++)) & 1;*//*#3*//*
    i = i << 1;
    line->data[24] = ((new_imm & i) >> (shift++)) & 1;*//*#4*//*
    i = i << 1;
    line->data[25] = ((new_imm & i) >> (shift++)) & 1;*//*#5*//*
    i = i << 1;
    line->data[26] = ((new_imm & i) >> (shift++)) & 1;*//*#6*//*
    i = i << 1;
    line->data[27] = ((new_imm & i) >> (shift++)) & 1;*//*#7*//*
    i = i << 1;
    line->data[28] = ((new_imm & i) >> (shift++)) & 1;*//*#8*//*
    i = i << 1;
    line->data[29] = ((new_imm & i) >> (shift++)) & 1;*//*#9*//*
    i = i << 1;
    line->data[30] = ((new_imm & i) >> (shift++)) & 1;*//*#10*//*
    i = i << 1;
    line->data[20] = ((new_imm & i) >> (shift++)) & 1;*//*#11*//*
    i = i << 1;
    line->data[12] = ((new_imm & i) >> (shift++)) & 1;*//*#12*//*
    i = i << 1;
    line->data[13] = ((new_imm & i) >> (shift++)) & 1;*//*#13*//*
    i = i << 1;
    line->data[14] = ((new_imm & i) >> (shift++)) & 1;*//*#14*//*
    i = i << 1;
    line->data[15] = ((new_imm & i) >> (shift++)) & 1;*//*#15*//*
    i = i << 1;
    line->data[16] = ((new_imm & i) >> (shift++)) & 1;*//*#16*//*
    i = i << 1;
    line->data[17] = ((new_imm & i) >> (shift++)) & 1;*//*#17*//*
    i = i << 1;
    line->data[18] = ((new_imm & i) >> (shift++)) & 1;*//*#18*//*
    i = i << 1;
    line->data[19] = ((new_imm & i) >> (shift++)) & 1;*//*#19*//*
    i = i << 1;
    line->data[31] = ((new_imm & (1 << 31)) >> 31) & 1;*//*#20*//*
}*/
