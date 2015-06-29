#include<stdio.h>
  /********************************************************************
   SUBROUTINE: read5_
   PURPOSE:  To read unformatted data from standard input
   ARGUMENTS:  ptrarray - a pointer to the record( in the form of an
                          array ), which will read len number of
                          unformatted data, len=number pointed to
                          by ptrlen.
   RETURNS : if no error during read : the number of records read
             if error occured during read : -1

   ********************************************************************/

        int read5_(ptrarray, ptrlen)
        int *ptrarray;
        int *ptrlen;
        {
          int irec;
          irec = fread(ptrarray,8,*ptrlen,stdin);
          if(irec==*ptrlen) {
             return(irec);
          }
          else {
             return(-1);
          }

        }




  /********************************************************************
   SUBROUTINE: write6_
   PURPOSE:  To write unformatted data to standard output
   ARGUMENTS:  ptrarray - a pointer to the record( in the form of an
                          array ), which contains the len number of
                          unformatted data points to be written out
                          , len=number pointed to by ptrlen

   RETURNS : the number of records written out

   ********************************************************************/

int write6_(ptrarray, ptrlen)
  int *ptrarray, *ptrlen;

{
  fwrite(ptrarray,8,*ptrlen,stdout);
  return(*ptrlen);
}


void vx2iei2 (buffer, n)
/* converts VAX I*2 to IEEE I*2, and vice-versa */
  unsigned char *buffer;
  int  *n;       /* no of I*2 to convert */
{

  int i;
  unsigned char b0, *temp;

  temp = buffer;

  for (i=0; i<*n; i++)
  {  b0       = temp [1];
     temp [1] = temp [0];
     temp [0] = b0;

     temp +=2;
  }

}



