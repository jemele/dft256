/*
This is DFT computation using matrix vector multiplication form.
INPUT:
	In_R, In_I[]: Real and Imag parts of Complex signal in time domain.
OUTPUT:
	In_R, In_I[]: Real and Imag parts of Complex signal in frequency domain.

*/

#include<stdio.h>
#include <stdlib.h>
#include<iostream>
#include <math.h>
#include "dft.h"

DTYPE In_R[SIZE], In_I[SIZE];

int main()
{
	FILE *fp;
	DTYPE pi = 3.141592653589;
	printf("INPUTS\n");
	for(int i=0; i<SIZE; i++)
	{
		In_R[i] = i;
		In_I[i] = 0.0;

	}
	


	//Perform DFT
	dft(In_R, In_I);
	//Print output
	fp=fopen("out.dat", "w");
	printf("\nPrinting DFT Output\n");
	for(int i=0; i<SIZE; i++){
		printf("%4d\t%f\t%f\n",i,In_R[i],In_I[i]);
		fprintf(fp, "%4d\t%f\t%f\n",i,In_R[i],In_I[i]);
	}
	fclose(fp);

	//Check against golden output.

  printf ("Comparing against output data \n");
  if (system("diff -w out.dat out.gold.dat")) {
	fprintf(stdout, "*******************************************\n");
	fprintf(stdout, "FAIL: Output DOES NOT match the golden output\n");
	fprintf(stdout, "*******************************************\n");
     return 1;
  } else {
	fprintf(stdout, "*******************************************\n");
	fprintf(stdout, "PASS: The output matches the golden output!\n");
	fprintf(stdout, "*******************************************\n");
     return 0;
  }

}
