/*
** simple error demonstration to demonstrate power of valgrind
** Julian M. Kunkel - 17.04.2008
*/

#include <stdio.h>
#include <stdlib.h>

int *
mistake1 ()
{
  int* buf = (int[6]){ 1, 1, 2, 3, 4, 5 };
  return buf;
}

int *
mistake2 ()
{
  int *buf = (int *) malloc (sizeof (int) * 3);
  buf[1] = 2;
  return buf;
}

int *
mistake3 ()
{
  //int mistake2_ = 0; 
  //int *buf = (int *) &mistake2;
  int *buf = (int *) malloc (sizeof (int));  
  buf[0] = 3;
  return buf;
}

int *
mistake4 ()
{
  int *buf = (int *) malloc (sizeof (int) * 5);
  buf[0] = 4;
  //free (buf);
  return buf;
}

int
main (void)
{
  /* Modifizieren Sie die folgende Zeile nicht */
  int *p[4] = { &mistake1 ()[1], &mistake2 ()[1], mistake3 (), mistake4 () };

  printf ("1 %d\n", *p[0]);
  printf ("2 %d\n", *p[1]);
  printf ("3 %d\n", *p[2]);
  printf ("4 %d\n", *p[3]);

  /* mhh muss hier noch etwas gefreed werden? */
  /* Fügen sie hier die korrekten aufrufe von free() ein */
  free (p[0]);
  free (p[1]);
  free (p[2]);
  free (p[3]);

  return 0;
}
