#ifndef COMPRESSION_H
#define COMPRESSION_H

/* Your code here... */
enum instruction_type {
    add,    /*add*/
    jalr,   /*jalr*/
    addi,   /*addi*/
    lui,    /*lui*/
    slli,   /*slli*/
    lw,     /*lw*/
    sw,     /*sw*/
    and,    /*and*/
    or,     /*or*/
    xor,    /*xor*/
    sub,    /*sub*/
    beq,    /*beq*/
    bne,    /*bne*/
    srli,   /*srli*/
    srai,   /*srai*/
    andi,   /*andi*/
    jal,    /*jal*/
    none,   /*none*/
    none_blt,   /*none_blt*/
    none_bge,   /*none_bge*/
    none_bltu,  /*none_bltu*/
    none_bgeu,   /*none_bgeu*/
    none_jal,    /*none_jar*/
    none_beq,    /*none_beq*/
    none_bne     /*none_bne*/
};

typedef struct {
    int data[32];   /*data*/
    int opcode;     /*opcode*/
    int rd;         /*rd*/
    int rs1;        /*rs1*/
    int rs2;        /*rs2*/
    int imm;        /*imm*/
    int shamt;      /*shamt*/
    int funct3;     /*funct3*/
    int funct7;     /*funct7*/
    enum instruction_type type;/*type*/
} instruction;

typedef struct {
    int data[16];   /*output data*/
    enum instruction_type type; /*type*/
} compressed_instruction;

/*check the type of instruction*/
void check_type(instruction *line);

/*check constrain of add*/
void check_add(instruction *line, int rd, int rs1, int rs2);

/*check constrain of lui*/
void check_lui(instruction *line, int rd, int imm);

/*check constrain of slli*/
void check_slli(instruction *line, int rd, int rs1);

/*check constrain of addi*/
void check_addi(instruction *line, int rd, int rs1, int imm);

/*check constrain of srli*/
void check_srli(instruction *line, int rd, int rs1);

/*check constrain of srai*/
void check_srai(instruction *line, int rd, int rs1);

/*check constrain of jalr*/
void check_jalr(instruction *line, int rd, int rs1, int imm);

/*check constrain of and*/
void check_and(instruction *line, int rd, int rs1, int rs2);

/*check constrain of or*/
void check_or(instruction *line, int rd, int rs1, int rs2);

/*check constrain of xor*/
void check_xor(instruction *line, int rd, int rs1, int rs2);

/*check constrain of sub*/
void check_sub(instruction *line, int rd, int rs1, int rs2);

/*check constrain of andi*/
void check_andi(instruction *line, int rd, int rs1, int imm);

/*check constrain of jal*/
void check_jal(instruction *line, int rd);

/*check constrain of lw*/
void check_lw(instruction *line, int rd, int rs1, int imm);

/*check constrain of sw*/
void check_sw(instruction *line, int rs1, int imm);

/*check constrain of beq*/
void check_beq(instruction *line, int rs1, int rs2);

/*check constrain of bne*/
void check_bne(instruction *line, int rs1, int rs2);

/*compress CR*/
compressed_instruction *compression_CR(instruction *line, compressed_instruction *_line);

/*compress CI*/
compressed_instruction *compression_CI(instruction *line, compressed_instruction *_line);

/*compress CL*/
compressed_instruction *compression_CL(instruction *line, compressed_instruction *_line);

/*compress CS*/
compressed_instruction *compression_CS(instruction *line, compressed_instruction *_line);

/*compress CB*/
compressed_instruction *
compression_CB(instruction *line, const int *address, int _line_No, compressed_instruction *_line);

/*compress CJ*/
compressed_instruction *
compression_CJ(instruction *line, const int *address, int _line_No, compressed_instruction *_line);

/*change offsets of bltu/bgeu/bge/blt*/
void change_offset(instruction *line, const int *address, int _line_No);

void change_offset_jal(instruction *line, const int *address, int _line_No);

#endif
