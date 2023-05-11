/*
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
Research Paper: Successful Simple Side Channel Analysis: 
Vulnerability of an atomic pattern kP algorithm implemented with a constant time crypto library to simple electromagnetic analysis attacks (accepted in MECO 2023)
  
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <flecc_in_c/types.h>
#include "flecc_in_c/io/io.h"
#include "flecc_in_c/bi/bi.h"
#include "flecc_in_c/eccp/eccp.h"
#include "flecc_in_c/utils/param.h"
#include "flecc_in_c/utils/parse.h"
#include "flecc_in_c/gfp/gfp.h"


//curve P-256
const char *curve_type_str = "secp256r1";
//mod p 
const char *Rsq_str = "4fffffffdfffffffffffffffefffffffbffffffff0000000000000003";

//coordinates
const char *q0x_str = "6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296";
const char *q0y_str = "4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5";
const char *q0z_str = "1";
const char *q0z1_str = "0";
const char *q0z2_str = "0";

//k scalar value
const char *KB =  "1000111011110111111101"; //23BDFD


const gfp_prime_data_t * dataWords;
eccp_parameters_t curve_params;
gfp_prime_data_t prime_data;

//Registers Initialization
gfp_t R0, R1, R2, R3;
gfp_t tempR0;

gfp_t Q0z;
gfp_t Q0z1;
gfp_t Q0z2;
gfp_t Q0x;
gfp_t Q0y;

gfp_t Px;
gfp_t Py;

gfp_t r_sq;

//Transition to Affinates Cordinates
gfp_t X_A;
gfp_t Y_A;
gfp_t Q0z1_temp;


void parse_bigint( const char *string, uint_t *big_int, const int bi_length ) {
    int len = strlen( string );
    bigint_parse_hex_var( big_int, bi_length, string, len );
}


void PointDoubling(gfp_t X1, gfp_t X2, gfp_t X3, gfp_t Z1, gfp_t Z2, gfp_prime_data_t p) {

    //OP 1   R0 <- X3 * X3
     gfp_mult_two_mont(R0, X3, X3, curve_params.prime_data.prime, r_sq);

    //OP 2   R2 <- X2 + X2
    gfp_gen_add(R2, X2, X2,curve_params.prime_data.prime ) ;
    
    //OP 3   R1 <- X1 - R0
    gfp_gen_subtract( R1, X1, R0, curve_params.prime_data.prime );

    //OP 4   Z1 <- X2 * R2
    gfp_mult_two_mont(Z1, X2, R2, curve_params.prime_data.prime,r_sq);

    //OP 5   X2 <- Z1 + Z1
    gfp_gen_add(X2, Z1, Z1, curve_params.prime_data.prime ) ;

    // dummy instructions to match the Point Addition//
    bigint_copy_var(tempR0, R0, curve_params.order_n_data.words );

    //OP 6   R3 <- R2 * X3
    gfp_mult_two_mont(R3, R2, X3, curve_params.prime_data.prime,r_sq);
    
    //OP 7   R2 <- X2 * X1
    gfp_mult_two_mont(R2, X2, X1, curve_params.prime_data.prime,r_sq);
    
    //OP 8   X1 <- X1 + R0
    // dummy instructions to match the Point Addition//
    gfp_gen_add(X1, X1, R0, curve_params.prime_data.prime ) ;

    //OP 9   R0 <- R1 * X1
    gfp_mult_two_mont(R0, R1, X1, curve_params.prime_data.prime,r_sq);

    //OP 10  R1 <- Z1 * X2
    gfp_mult_two_mont(R1, Z1, X2, curve_params.prime_data.prime,r_sq);
    
    //OP 11  X1 <- R0 + R0
    gfp_gen_add(X1, R0, R0, curve_params.prime_data.prime ) ;
    
    //OP 12  R0 <- R0 + X1
    gfp_gen_add(R0, R0, X1, curve_params.prime_data.prime ) ;
 
    //OP 13  X1 <- R0^2
    gfp_mult_two_mont(X1, R0, R0, curve_params.prime_data.prime,r_sq);

    //OP 14  X1 <- X1 - R2
    gfp_gen_subtract( X1, X1, R2, curve_params.prime_data.prime );

    //OP 15  Z1 <- R3^2
    gfp_mult_two_mont(Z1, R3, R3, curve_params.prime_data.prime,r_sq);

    //OP 16  X1 <- X1 - R2
    gfp_gen_subtract( X1, X1, R2, curve_params.prime_data.prime );
 
    //OP 17  R2 <- R2 - X1
    gfp_gen_subtract( R2, R2, X1, curve_params.prime_data.prime );
 
    //OP 18  Z2 <- Z1 * R3
    gfp_mult_two_mont(Z2, Z1, R3, curve_params.prime_data.prime,r_sq);

    //OP 19  X2 <- R0 * R2
    gfp_mult_two_mont(X2, R0, R2, curve_params.prime_data.prime,r_sq);
 
    //OP 20  X3 <- R3
    bigint_copy_var(X3, R3, curve_params.order_n_data.words);
 
    //OP 21  X2 <- X2 - R1
    gfp_gen_subtract( X2, X2, R1, curve_params.prime_data.prime );
 
}

void PointAddition(gfp_t X1, gfp_t X2, gfp_t X3, gfp_t Z1, gfp_t Z2, gfp_prime_data_t p, gfp_t X, gfp_t Y) {

    //OP 1   R1 <- X * Z1
    gfp_mult_two_mont(R1, X, Z1, curve_params.prime_data.prime,r_sq);

    //OP 2   R2 <- X2 + X2
    gfp_gen_add( R2, X2, X2, curve_params.prime_data.prime ) ;

    //OP 3   R1 <- R1 - X1
    gfp_gen_subtract( R1, R1, X1, curve_params.prime_data.prime ) ;
  
    //OP 4   R2 <- R1 * R1
    gfp_mult_two_mont(R2, R1, R1, curve_params.prime_data.prime,r_sq);
    
    //OP 5   R0 <- R2 + R2
    gfp_gen_add( R0, R2, R2, curve_params.prime_data.prime ) ;
    bigint_copy_var(tempR0, R0, curve_params.order_n_data.words );

    //OP 6   R3 <- X1 * R2
    gfp_mult_two_mont(R3, X1, R2, curve_params.prime_data.prime,r_sq);
      
    //OP 7   R0 <- Y * Z2
    gfp_mult_two_mont(R0, Y, Z2, curve_params.prime_data.prime,r_sq);
   
    //OP 8   Z2 <- Z2 + R0
    gfp_gen_add( Z2, Z2, tempR0, curve_params.prime_data.prime );
       
    //OP 9   Z2 <- R1 * R2
    gfp_mult_two_mont(Z2, R1, R2, curve_params.prime_data.prime,r_sq);
   
    //OP 10  R2 <- X3 * R1
    gfp_mult_two_mont(R2, X3, R1, curve_params.prime_data.prime,r_sq);

    //OP 11  X1 <- R3 + R3
    gfp_gen_add(X1, R3, R3, curve_params.prime_data.prime ) ;
    
    //OP 12  X1 <- Z2 + X1
    gfp_gen_add(X1, Z2, X1, curve_params.prime_data.prime ) ;
       
    //OP 13  Z1 <- X1^2
    gfp_mult_two_mont(Z1, X1, X1, curve_params.prime_data.prime,r_sq);
    
    //OP 14  R0 <- R0 - X2
    gfp_gen_subtract( R0, R0, X2, curve_params.prime_data.prime );
    
    //OP 15  R1 <- R0^2
    gfp_mult_two_mont(R1, R0, R0, curve_params.prime_data.prime,r_sq);
    
    //OP 16  X1 <- R1 - X1
    gfp_gen_subtract( X1, R1, X1, curve_params.prime_data.prime );
    
    //OP 17  R1 <- R3 - X1
    gfp_gen_subtract( R1, R3, X1, curve_params.prime_data.prime );
    
    //OP 18  R3 <- R1 * R0
    gfp_mult_two_mont(R3, R1, R0, curve_params.prime_data.prime,r_sq);
    
    //OP 19  R0 <- X2 * Z2
    gfp_mult_two_mont(R0, X2, Z2, curve_params.prime_data.prime,r_sq);
    
    //OP 20  X3 <- R2
    bigint_copy_var( X3, R2, curve_params.order_n_data.words );
    
    //OP 21  X2 <- R3 - R0
    gfp_gen_subtract( X2, R3, R0, curve_params.prime_data.prime );
   
}


int main(void)
{

       // load curve
       curve_params.curve_type = param_get_curve_type_from_name( curve_type_str, strlen( curve_type_str ) );
       param_load( &curve_params, curve_params.curve_type );

       prime_data = curve_params.order_n_data;
       dataWords = curve_params.order_n_data.words;


       parse_bigint( q0z_str, Q0z, curve_params.order_n_data.words);
       parse_bigint( q0z1_str, Q0z1, curve_params.order_n_data.words);
       parse_bigint( q0z2_str, Q0z2, curve_params.order_n_data.words);
       parse_bigint( q0x_str, Q0x, curve_params.order_n_data.words);
       parse_bigint( q0y_str, Q0y, curve_params.order_n_data.words);
       parse_bigint( q0x_str, Px, curve_params.order_n_data.words);
       parse_bigint( q0y_str, Py, curve_params.order_n_data.words);

       parse_bigint( Rsq_str, r_sq, curve_params.order_n_data.words);


       int l = strlen(KB);

       //Step1 - Q = P
       bigint_copy_var( Px, Q0x, curve_params.order_n_data.words );
       
       //Py = Q0y;
       bigint_copy_var( Py, Q0y, curve_params.order_n_data.words);
       
       //Step2
       int t = l - 1;
       int i = t;
       
       for ( i = t; i > 0; i--) {
          
           //Step3
           PointDoubling(Q0x, Q0y, Q0z, Q0z1, Q0z2, prime_data);
          
           if (KB[l - i] == '1') {
               PointAddition(Q0x, Q0y, Q0z, Q0z1, Q0z2, prime_data, Px, Py);          
           }
           else if (KB[l - i] != '0') {              
               break;
           }

       }
	   
	   //Transition to Affinates Coordinates
	  if (KB[l - 1] == '1') {
		  
		   //X3^2
           gfp_mult_two_mont(Q0z1, Q0z, Q0z, curve_params.prime_data.prime,r_sq);
           bigint_copy_var(Q0z1_temp, Q0z1, curve_params.order_n_data.words );	
		   //X_A = 1/X3^2
           gfp_binary_euclidean_inverse( X_A, Q0z1_temp, curve_params.prime_data.prime);   
		   //X_A = X1/X3^2
           gfp_mult_two_mont(X_A, Q0x, X_A, curve_params.prime_data.prime,r_sq);
		   
		   //X3^3
           gfp_mult_two_mont(Q0z2, Q0z, Q0z1, curve_params.prime_data.prime,r_sq);
		   //Y_A = 1/X3^3
           gfp_binary_euclidean_inverse( Y_A, Q0z2, curve_params.prime_data.prime );
		   //Y_A = X2/X3^2
           gfp_mult_two_mont(Y_A, Q0y, Y_A, curve_params.prime_data.prime,r_sq);
           
		 
           printf("X_A1: ");
           io_print_bigint_var(X_A, curve_params.order_n_data.words);
           printf("Y_A1: ");
           io_print_bigint_var( Y_A, curve_params.order_n_data.words);
        }
        else if (KB[l - 1] == '0') {

		   //X_A = 1/Z1
           gfp_binary_euclidean_inverse( X_A, Q0z1, curve_params.prime_data.prime);
		   //X_A = X1/Z1
           gfp_mult_two_mont(X_A, Q0x, Q0z1, curve_params.prime_data.prime,r_sq);

		   //Y_A = 1/Z2
           gfp_binary_euclidean_inverse( Y_A, Q0z2, curve_params.prime_data.prime);
		   //Y_A = X2/Z2
           gfp_mult_two_mont(Y_A, Q0y, Q0z2, curve_params.prime_data.prime,r_sq);

           printf("X_A0: ");
           io_print_bigint_var(X_A, curve_params.order_n_data.words);
           printf("Y_A0: ");
           io_print_bigint_var( Y_A, curve_params.order_n_data.words);
        }

        else printf("Error!\n");


    return 0;
}
