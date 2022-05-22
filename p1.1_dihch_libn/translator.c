/*  Project 1.1: RISC-V instructions to RISC-V compressed instructions in C89.
    The following is the starter code provided for you. To finish the task, you 
    should define and implement your own functions in translator.c, compression.c, 
    utils.c and their header files.
    Please read the problem description before you start.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "src/compression.h"
#include "src/utils.h"

#include "translator.h"

/*check if file can be correctly opened */
static int open_files(FILE **input, FILE **output, const char *input_name, const char *output_name) {
    *input = fopen(input_name, "r");
    if (!*input) { /* open input file failed */
        printf("Error: unable to open input file: %s\n", input_name);
        return -1;
    }
    *output = fopen(output_name, "w");
    if (!*output) { /* open output file failed */
        printf("Error: unable to open output file: %s\n", output_name);
        fclose(*input);
        return -1;
    }
    return 0;
    /* no problem opening files */
}

static int close_files(FILE **input, FILE **output) {
    fclose(*input);
    fclose(*output); /* close the files at the end */
    return 0;
}

static void print_usage_and_exit() {
    printf("Usage:\n");
    printf("Run program with translator <input file> <output file>\n"); /* print the correct usage of the program */
    exit(0);
}


/*Run the translator 
*/
int translate(const char *in, const char *out) {
    FILE *input, *output;

    int i, j, status = 1, compress_address[51];
    /*i,j iterators*/

    compressed_instruction *res = (compressed_instruction *) malloc(sizeof(compressed_instruction));
    /*result of compression*/

    instruction lines[51];
    /*all input lines*/

    for (i = 0; i < 51; ++i) {
        lines[i].data[0] = -1;
    }

    memset(res->data,0,sizeof(compressed_instruction));
    /*initialization of res->data*/


    memset(compress_address, 0, 51 * sizeof(int));
    /*initialization of lines*/

    /*edge case*/
    if (!in || open_files(&input, &output, in, out) != 0) {
        free(res);
        return 1;
    }

    /*read data and preprocess*/
    for (i = 0; i < 51 && status; ++i) {
        status = read_line(input, lines[i].data);
        check_type(&lines[i]);
        /*classify the type of each line*/
        if (lines[i].type == none || lines[i].type == none_blt || lines[i].type == none_bge ||
            lines[i].type == none_bgeu || lines[i].type == none_bltu || lines[i].type == none_jal ||
            lines[i].type == none_beq || lines[i].type == none_bne) {
            /*if the line can not be compressed*/
            compress_address[i + 1] = compress_address[i] + 2;
            /*use 16-bit as a unit, so we add 2 unit to represent offset = 32-bit*/
        } else compress_address[i + 1] = compress_address[i] + 1;
        /*use 16-bit as a unit, so we plus 1 unit*/
    }

    j = i;
    /*remember the number of lines*/

    for (i = 0; i < j && lines[i].data[0] != -1 ; ++i) {
        enum instruction_type type = lines[i].type;
        if (type == add || type == jalr) {                                                          /*CR*/
            res = compression_CR(&lines[i], res);
            write_line(output, res->data, 1);
        } else if (type == addi || type == lui || type == slli) {                                   /*CI*/
            res = compression_CI(&lines[i], res);
            write_line(output, res->data, 1);
        } else if (type == lw) { /*CL*/
            res = compression_CL(&lines[i], res);
            write_line(output, res->data, 1);
        } else if (type == sw || type == and || type == or || type == xor || type == sub) {         /*CS*/
            res = compression_CS(&lines[i], res);
            write_line(output, res->data, 1);
        } else if (type == beq || type == bne || type == srli || type == srai || type == andi) {    /*CB*/
            res = compression_CB(&lines[i], compress_address,i, res);
            write_line(output, res->data, 1);
        } else if (type == jal) {                                                                   /*CJ*/
            res = compression_CJ(&lines[i], compress_address,i, res);
            write_line(output, res->data, 1);
        } else if (type == none_bltu || type == none_bgeu ||
                   type == none_bge || type == none_blt ||
                   type == none_beq || type == none_bne) {          /*change offsets in bltu/bgeu/bge/blt*/
            change_offset(&lines[i], compress_address,i);
            write_line(output, lines[i].data, 0);
        }
        else if (type == none_jal) write_line(output,lines[i].data,0);
        /*else if (type == none_jal) {
            none_jal
            change_offset_jal(&lines[i], compress_address, i);
            write_line(output, lines[i].data, 0);
        }*/
        else if (type == none) { /*none*/                                                          /*none*/
            write_line(output, lines[i].data, 0);
        }
    }
    /*free the allocated memory*/
    free(res);

    close_files(&input, &output);
    /*close I/O files*/
    return 0;
}

/* main func */
int main(int argc, char **argv) {
    char *input_fname, *output_fname;
    /*I/O*/
    int err;

    if (argc != 3)
        /* need correct arguments */
        print_usage_and_exit();

    input_fname = argv[1]; /*argument 1*/
    output_fname = argv[2];/*argument 2*/

    err = translate(input_fname, output_fname); /* main translation process */
    if (err)
        printf("One or more errors encountered during translation operation.\n"); /* something wrong */
    else
        printf("Translation process completed successfully.\n"); /* correctly output */

    return 0;
}
