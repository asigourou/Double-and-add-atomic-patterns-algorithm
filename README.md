# Double-and-add-atomic-patterns-algorithm
Our software, is an implementation of a kP operation for Elliptic Curves over GF(p) based on the atomicity principle of Rondepierre. 
The source code is written in the programming language C using the constant run time functions of the open-source library FLECC_IN_C .
This version of the software loads first the secp256r1 curve and takes specific coordinates of the curve and a 22-bit key as example inputs.
In the main function implements the "left-to-right double and add algorithm" based on the atomicity principle,
where Point Doubling(PD) and Point Addition(PA) functions are constructed by the same sequence of instructions.
The PD and PA functions are created based on the steps of the algorithm published in Table 2 of the paper:
"Ievgen Kabin, Zoya Dyka and Peter Langendoerfer: Atomicity and Regularity Principles Do Not Ensure Full Resistance of ECC Designs against Single-Trace Attacks".
This release is intended for developers and researchers interested in implementing kP operations for Elliptic Curves over GF(p).


Authors: Alkistis Aikaterini Sigourou, Ievgen Kabin 
{sigourou, kabin}@ihp-microelectronics.com

IHP—Leibniz-Institut für Innovative Mikroelektronik, 15236 Frankfurt, Germany

Open Source Library: FLECC_IN_C: https://github.com/IAIK/flecc_in_c

This software has been created for reasearch purposes, in the limits of testing the resistance of the algorithm on Side Channel Attacks.

To meet our purpose we have modified the function gfp_mult_two_mont.

/////////////////////////////////////////////  Original Version*   /////////////////////////////////////////

void gfp_mult_two_mont( gfp_t res, const gfp_t a, const gfp_t b, const gfp_prime_data_t *prime_data,gfp_t r_squared) {

    gfp_mont_multiply( res, a, b, prime_data );
    
    gfp_mont_multiply( res, res, prime_data->r_squared, prime_data );
}

*https://github.com/IAIK/flecc_in_c/blob/develop/src/gfp/gfp_mont.c (last accessed on 27/04/2023)


/////////////////////////////////////////////  Our Version  ///////////////////////////////////////////////////

void gfp_mult_two_mont( gfp_t res, const gfp_t a, const gfp_t b, const gfp_prime_data_t *prime_data,gfp_t r_squared) {

    gfp_mont_multiply( res, a, b, prime_data );
    
    gfp_mont_multiply( res, res, r_squared, prime_data );
}

file : url

