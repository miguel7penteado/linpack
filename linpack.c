/*
**
** LINPACK.C        Linpack benchmark, calculates FLOPS.
**                  (FLoating Point Operations Per Second)
**
** Translated to C by Bonnie Toy 5/88
**
** Modified by Will Menninger, 10/93, with these features:
**  (modified on 2/25/94  to fix a problem with daxpy  for
**   unequal increments or equal increments not equal to 1.
**     Jack Dongarra)
**
** - Defaults to double precision.
** - Averages ROLLed and UNROLLed performance.
** - User selectable array sizes.
** - Automatically does enough repetitions to take at least 10 CPU seconds.
** - Prints machine precision.
** - ANSI prototyping.
**
** To compile:  cc -O -o linpack linpack.c -lm
**
**
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>

#define DP

#ifdef SP
#define ZERO        0.0
#define ONE         1.0
#define PREC        "Single"
#define BASE10DIG   FLT_DIG

typedef float   TIPO_REAL;
#endif

#ifdef DP
#define ZERO        0.0e0
#define ONE         1.0e0
#define PREC        "Double"
#define BASE10DIG   DBL_DIG

typedef double  TIPO_REAL;
#endif

// declaração das funções 

static TIPO_REAL linpack  (long nreps,int arsize);
static void      matgen   (TIPO_REAL *a,int lda,int n,TIPO_REAL *b,TIPO_REAL *norma);
static void      dgefa    (TIPO_REAL *a,int lda,int n,int *ipvt,int *info,int roll);
static void      dgesl    (TIPO_REAL *a,int lda,int n,int *ipvt,TIPO_REAL *b,int job,int roll);
static void      daxpy_r  (int n,TIPO_REAL da,TIPO_REAL *dx,int incx,TIPO_REAL *dy,int incy);
static TIPO_REAL ddot_r   (int n,TIPO_REAL *dx,int incx,TIPO_REAL *dy,int incy);
static void      dscal_r  (int n,TIPO_REAL da,TIPO_REAL *dx,int incx);
static void      daxpy_ur (int n,TIPO_REAL da,TIPO_REAL *dx,int incx,TIPO_REAL *dy,int incy);
static TIPO_REAL ddot_ur  (int n,TIPO_REAL *dx,int incx,TIPO_REAL *dy,int incy);
static void      dscal_ur (int n,TIPO_REAL da,TIPO_REAL *dx,int incx);
static int       idamax   (int n,TIPO_REAL *dx,int incx);
static TIPO_REAL second   (void);

static void *mempool;


int main(void)
{
    char    *entrada_tamanhoMatriz;
    int     arsize;
    long    arsize2d,nreps;
    size_t  malloc_arg,memreq;

    entrada_tamanhoMatriz = getenv("LINPACK_TAMANHO_MATRIZ");
	
    if (entrada_tamanhoMatriz == NULL) 
    {
        arsize = 200;
    } else 
    {
        arsize = atoi(entrada_tamanhoMatriz);
    }

	arsize/=2;
	arsize*=2;
	
	if (arsize<10)
	{
		printf("Muito pequena.\n");
		return 1;
	}

	arsize2d = (long)arsize*(long)arsize;
	memreq=arsize2d*sizeof(TIPO_REAL)+(long)arsize*sizeof(TIPO_REAL)+(long)arsize*sizeof(int);
	printf("Memoria requisitada:  %ldK.\n",(memreq+512L)>>10);
	malloc_arg=(size_t)memreq;

	if (malloc_arg!=memreq || (mempool=malloc(malloc_arg))==NULL)
	{
		printf("A memoria é insuficiente para o tamanho da matriz fornecido.\n\n");
		return 2;
	}
	printf("\n\nLINPACK benchmark, %s precision.\n",PREC);
	printf("Precisao da maquina:  %d digits.\n",BASE10DIG);
	printf("Matriz de tamanho %d X %d.\n",arsize,arsize);
	printf("Desempenho médio enrolado e desenrolado:\n\n");
	printf("    Reps Time(s) DGEFA   DGESL  OVERHEAD    KFLOPS\n");
	printf("----------------------------------------------------\n");
	nreps=1;

	while (linpack(nreps,arsize)<10.)
	{
		nreps*=2;
	}

	free(mempool);
	printf("\n");
	
}


static TIPO_REAL linpack(long nreps,int arsize)
{
	TIPO_REAL  *a,*b;
	TIPO_REAL   norma,t1,kflops,tdgesl,tdgefa,totalt,toverhead,ops;
	int   *ipvt,n,info,lda;
	long   i,arsize2d;
	
	lda = arsize;
	n = arsize/2;
	arsize2d = (long)arsize*(long)arsize;
	ops=((2.0*n*n*n)/3.0+2.0*n*n);
	a=(TIPO_REAL *)mempool;
	b=a+arsize2d;
	ipvt=(int *)&b[arsize];
	tdgesl=0;
	tdgefa=0;
	totalt=second();

    for (i=0;i<nreps;i++)
	{
		matgen(a,lda,n,b,&norma);
		t1 = second();
		dgefa(a,lda,n,ipvt,&info,1);
		tdgefa += second()-t1;
		t1 = second();
		dgesl(a,lda,n,ipvt,b,0,1);
		tdgesl += second()-t1;
	}

    for (i=0;i<nreps;i++)
	{
		matgen(a,lda,n,b,&norma);
		t1 = second();
		dgefa(a,lda,n,ipvt,&info,0);
		tdgefa += second()-t1;
		t1 = second();
		dgesl(a,lda,n,ipvt,b,0,0);
		tdgesl += second()-t1;
	}
    totalt=second()-totalt;

    if (totalt<0.5 || tdgefa+tdgesl<0.2)
	{
        return(0.);
	}
    kflops=2.*nreps*ops/(1000.*(tdgefa+tdgesl));
    toverhead=totalt-tdgefa-tdgesl;

    if (tdgefa<0.)
	{
        tdgefa=0.;
	}

    if (tdgesl<0.)
	{
        tdgesl=0.;
	}

    if (toverhead<0.)
	{
        toverhead=0.;
	}

    printf("%8ld %6.2f %6.2f%% %6.2f%% %6.2f%%  %9.3f\n",
            nreps,totalt,100.*tdgefa/totalt,
            100.*tdgesl/totalt,100.*toverhead/totalt,
            kflops);

    return(totalt);
}


/*
** For matgen,
** We would like to declare a[][lda], but c does not allow it.  In this
** function, references to a[i][j] are written a[lda*i+j].
*/
static void matgen(TIPO_REAL *a,int lda,int n,TIPO_REAL *b,TIPO_REAL *norma)

    {
    int init,i,j;

    init = 1325;
    *norma = 0.0;
    for (j = 0; j < n; j++)
        for (i = 0; i < n; i++)
            {
            init = (int)((long)3125*(long)init % 65536L);
            a[lda*j+i] = (init - 32768.0)/16384.0;
            *norma = (a[lda*j+i] > *norma) ? a[lda*j+i] : *norma;
            }
    for (i = 0; i < n; i++)
        b[i] = 0.0;
    for (j = 0; j < n; j++)
        for (i = 0; i < n; i++)
            b[i] = b[i] + a[lda*j+i];
    }


/*

DGEFA benchmark

We would like to declare a[][lda], but c does not allow it.  In this
function, references to a[i][j] are written a[lda*i+j].

  A função dgefa fatora uma matriz de precisão dupla por eliminação gaussiana.

  dgefa is usually called by dgeco, but it can be called
  directly with a saving in time if  rcond  is not needed.
  (time for dgeco) = (1 + 9/n)*(time for dgefa) .

  na entrada

     a       TIPO_REAL precision[n][lda]
             a matriz a ser fatorada.

     lda     inteiro
             a dimensão principal da matriz a.

     n       inteiro
             a ordem da matriz a.

  no retorno

     a       uma matriz triangular superior e os multiplicadores
             que foram usados para obtê-lo.
             a fatoração pode ser escrita  a = l*u  onde
             l  é um produto de permutação e unidade inferior
             matrizes triangulares e u é triangular superior.

     ipvt    integer[n]
             um vetor inteiro de índices de pivô.

     info    integer
             = 0  normal value.
             = k  if  u[k][k] .eq. 0.0 .  isso não é um erro
                  condição para esta sub-rotina, mas não
                  indica que dgesl ou dgedi serão divididos por zero
                  se chamado. use rcond em dgeco para um confiável
                  indicação de singularidade.

  linpack. this version dated 08/14/78 .
  cleve moler, Universidade do Novo México, laboratório nacional de Argonne.

  funções

  blas daxpy, dscal, idamax

*/
static void dgefa(TIPO_REAL *a,int lda,int n,int *ipvt,int *info,int roll)

    {
    TIPO_REAL t;
    int idamax(),j,k,kp1,l,nm1;

    /* eliminação gaussiana com pivotamento parcial */

    if (roll)
        {
        *info = 0;
        nm1 = n - 1;
        if (nm1 >=  0)
            for (k = 0; k < nm1; k++)
                {
                kp1 = k + 1;

                /* find l = pivot index */

                l = idamax(n-k,&a[lda*k+k],1) + k;
                ipvt[k] = l;

                /* pivô zero implica que esta coluna já está triangularizada */

                if (a[lda*k+l] != ZERO)
                    {

                    /* troque se necessário */

                    if (l != k)
                        {
                        t = a[lda*k+l];
                        a[lda*k+l] = a[lda*k+k];
                        a[lda*k+k] = t;
                        }

                    /* calcular multiplicadores */

                    t = -ONE/a[lda*k+k];
                    dscal_r(n-(k+1),t,&a[lda*k+k+1],1);

                    /* eliminação de linha com indexação de coluna */

                    for (j = kp1; j < n; j++)
                        {
                        t = a[lda*j+l];
                        if (l != k)
                            {
                            a[lda*j+l] = a[lda*j+k];
                            a[lda*j+k] = t;
                            }
                        daxpy_r(n-(k+1),t,&a[lda*k+k+1],1,&a[lda*j+k+1],1);
                        }
                    }
                else
                    (*info) = k;
                }
        ipvt[n-1] = n-1;
        if (a[lda*(n-1)+(n-1)] == ZERO)
            (*info) = n-1;
        }
    else
        {
        *info = 0;
        nm1 = n - 1;
        if (nm1 >=  0)
            for (k = 0; k < nm1; k++)
                {
                kp1 = k + 1;

                /* find l = pivot index */

                l = idamax(n-k,&a[lda*k+k],1) + k;
                ipvt[k] = l;

                /* pivô zero implica que esta coluna já está triangularizada */

                if (a[lda*k+l] != ZERO)
                    {

                    /* troque se necessário */

                    if (l != k)
                        {
                        t = a[lda*k+l];
                        a[lda*k+l] = a[lda*k+k];
                        a[lda*k+k] = t;
                        }

                    /* calcular multiplicadores */

                    t = -ONE/a[lda*k+k];
                    dscal_ur(n-(k+1),t,&a[lda*k+k+1],1);

                    /* eliminação de linha com indexação de coluna */

                    for (j = kp1; j < n; j++)
                        {
                        t = a[lda*j+l];
                        if (l != k)
                            {
                            a[lda*j+l] = a[lda*j+k];
                            a[lda*j+k] = t;
                            }
                        daxpy_ur(n-(k+1),t,&a[lda*k+k+1],1,&a[lda*j+k+1],1);
                        }
                    }
                else
                    (*info) = k;
                }
        ipvt[n-1] = n-1;
        if (a[lda*(n-1)+(n-1)] == ZERO)
            (*info) = n-1;
        }
    }


/*

DGESL benchmark

We would like to declare a[][lda], but c does not allow it.  In this
function, references to a[i][j] are written a[lda*i+j].

  dgesl resolve o sistema de precisão dupla
  a * x = b  or  trans(a) * x = b
  usando os fatores calculados por dgeco ou dgefa.

  na entrada

     a       double precision[n][lda]
             a saída de dgeco ou dgefa.

     lda     integer
             a dimensão principal da matriz a.

     n       integer
             the order of the matrix  a .

     ipvt    integer[n]
             o vetor pivô de dgeco ou dgefa.

     b       double precision[n]
             o vetor do lado direito.

     job     integer
             = 0         para resolver  a*x = b ,
             = nonzero   para resolver  trans(a)*x = b  onde
                         trans(a)  é a transposição.

 no retorno

     b       the solution vector  x .

  error condition

     a division by zero will occur if the input factor contains a
     zero on the diagonal.  technically this indicates singularity
     but it is often caused by improper arguments or improper
     setting of lda .  it will not occur if the subroutines are
     called correctly and if dgeco has set rcond .gt. 0.0
     or dgefa has set info .eq. 0 .

  to compute  inverse(a) * c  where  c  is a matrix
  with  p  columns
        dgeco(a,lda,n,ipvt,rcond,z)
        if (!rcond is Muito pequena){
             for (j=0,j<p,j++)
                     dgesl(a,lda,n,ipvt,c[j][0],0);
        }

  linpack. this version dated 08/14/78 .
  cleve moler, university of new mexico, argonne national lab.

  functions

  blas daxpy,ddot
*/
static void dgesl(TIPO_REAL *a,int lda,int n,int *ipvt,TIPO_REAL *b,int job,int roll)

    {
    TIPO_REAL    t;
    int     k,kb,l,nm1;

    if (roll)
        {
        nm1 = n - 1;
        if (job == 0)
            {

            /* job = 0 , solve  a * x = b   */
            /* first solve  l*y = b         */

            if (nm1 >= 1)
                for (k = 0; k < nm1; k++)
                    {
                    l = ipvt[k];
                    t = b[l];
                    if (l != k)
                        {
                        b[l] = b[k];
                        b[k] = t;
                        }
                    daxpy_r(n-(k+1),t,&a[lda*k+k+1],1,&b[k+1],1);
                    }

            /* now solve  u*x = y */

            for (kb = 0; kb < n; kb++)
                {
                k = n - (kb + 1);
                b[k] = b[k]/a[lda*k+k];
                t = -b[k];
                daxpy_r(k,t,&a[lda*k+0],1,&b[0],1);
                }
            }
        else
            {

            /* job = nonzero, solve  trans(a) * x = b  */
            /* first solve  trans(u)*y = b             */

            for (k = 0; k < n; k++)
                {
                t = ddot_r(k,&a[lda*k+0],1,&b[0],1);
                b[k] = (b[k] - t)/a[lda*k+k];
                }

            /* now solve trans(l)*x = y     */

            if (nm1 >= 1)
                for (kb = 1; kb < nm1; kb++)
                    {
                    k = n - (kb+1);
                    b[k] = b[k] + ddot_r(n-(k+1),&a[lda*k+k+1],1,&b[k+1],1);
                    l = ipvt[k];
                    if (l != k)
                        {
                        t = b[l];
                        b[l] = b[k];
                        b[k] = t;
                        }
                    }
            }
        }
    else
        {
        nm1 = n - 1;
        if (job == 0)
            {

            /* job = 0 , solve  a * x = b   */
            /* first solve  l*y = b         */

            if (nm1 >= 1)
                for (k = 0; k < nm1; k++)
                    {
                    l = ipvt[k];
                    t = b[l];
                    if (l != k)
                        {
                        b[l] = b[k];
                        b[k] = t;
                        }
                    daxpy_ur(n-(k+1),t,&a[lda*k+k+1],1,&b[k+1],1);
                    }

            /* now solve  u*x = y */

            for (kb = 0; kb < n; kb++)
                {
                k = n - (kb + 1);
                b[k] = b[k]/a[lda*k+k];
                t = -b[k];
                daxpy_ur(k,t,&a[lda*k+0],1,&b[0],1);
                }
            }
        else
            {

            /* job = nonzero, solve  trans(a) * x = b  */
            /* first solve  trans(u)*y = b             */

            for (k = 0; k < n; k++)
                {
                t = ddot_ur(k,&a[lda*k+0],1,&b[0],1);
                b[k] = (b[k] - t)/a[lda*k+k];
                }

            /* now solve trans(l)*x = y     */

            if (nm1 >= 1)
                for (kb = 1; kb < nm1; kb++)
                    {
                    k = n - (kb+1);
                    b[k] = b[k] + ddot_ur(n-(k+1),&a[lda*k+k+1],1,&b[k+1],1);
                    l = ipvt[k];
                    if (l != k)
                        {
                        t = b[l];
                        b[l] = b[k];
                        b[k] = t;
                        }
                    }
            }
        }
    }



/*
 Constante vezes um vetor mais um vetor.
 Jack Dongarra, linpack, 3/11/78.
 Versão ROLLED
*/
static void daxpy_r(int n,TIPO_REAL da,TIPO_REAL *dx,int incx,TIPO_REAL *dy,int incy)

    {
    int i,ix,iy;

    if (n <= 0)
        return;
    if (da == ZERO)
        return;

    if (incx != 1 || incy != 1)
        {

        /* código para incrementos desiguais ou incrementos iguais != 1 */

        ix = 1;
        iy = 1;
        if(incx < 0) ix = (-n+1)*incx + 1;
        if(incy < 0)iy = (-n+1)*incy + 1;
        for (i = 0;i < n; i++)
            {
            dy[iy] = dy[iy] + da*dx[ix];
            ix = ix + incx;
            iy = iy + incy;
            }
        return;
        }

    /* código para ambos os incrementos iguais a 1 */

    for (i = 0;i < n; i++)
        dy[i] = dy[i] + da*dx[i];
    }


/*
 Forma o produto escalar de dois vetores.
 Jack Dongarra, linpack, 3/11/78.
 Versão ROLLED
*/
static TIPO_REAL ddot_r(int n,TIPO_REAL *dx,int incx,TIPO_REAL *dy,int incy)

    {
    TIPO_REAL dtemp;
    int i,ix,iy;

    dtemp = ZERO;

    if (n <= 0)
        return(ZERO);

    if (incx != 1 || incy != 1)
        {

        /* code for unequal increments or equal increments != 1 */

        ix = 0;
        iy = 0;
        if (incx < 0) ix = (-n+1)*incx;
        if (incy < 0) iy = (-n+1)*incy;
        for (i = 0;i < n; i++)
            {
            dtemp = dtemp + dx[ix]*dy[iy];
            ix = ix + incx;
            iy = iy + incy;
            }
        return(dtemp);
        }

    /* código para ambos os incrementos iguais a 1 */

    for (i=0;i < n; i++)
        dtemp = dtemp + dx[i]*dy[i];
    return(dtemp);
    }


/*
 Escala um vetor por uma constante.
 Jack Dongarra, linpack, 3/11/78.
 Versão ROLLED
*/
static void dscal_r(int n,TIPO_REAL da,TIPO_REAL *dx,int incx)

    {
    int i,nincx;

    if (n <= 0)
        return;
    if (incx != 1)
        {

        /* código para incremento não igual a 1 */

        nincx = n*incx;
        for (i = 0; i < nincx; i = i + incx)
            dx[i] = da*dx[i];
        return;
        }

    /* código para incremento igual a 1 */

    for (i = 0; i < n; i++)
        dx[i] = da*dx[i];
    }


/*
** constante vezes um vetor mais um vetor.
** Jack Dongarra, linpack, 3/11/78.
** UNVersão ROLLED
*/
static void daxpy_ur(int n,TIPO_REAL da,TIPO_REAL *dx,int incx,TIPO_REAL *dy,int incy)

    {
    int i,ix,iy,m;

    if (n <= 0)
        return;
    if (da == ZERO)
        return;

    if (incx != 1 || incy != 1)
        {

        /* código para incrementos desiguais ou incrementos iguais != 1 */

        ix = 1;
        iy = 1;
        if(incx < 0) ix = (-n+1)*incx + 1;
        if(incy < 0)iy = (-n+1)*incy + 1;
        for (i = 0;i < n; i++)
            {
            dy[iy] = dy[iy] + da*dx[ix];
            ix = ix + incx;
            iy = iy + incy;
            }
        return;
        }

    /* código para ambos os incrementos iguais a 1 */

    m = n % 4;
    if ( m != 0)
        {
        for (i = 0; i < m; i++)
            dy[i] = dy[i] + da*dx[i];
        if (n < 4)
            return;
        }
    for (i = m; i < n; i = i + 4)
        {
        dy[i] = dy[i] + da*dx[i];
        dy[i+1] = dy[i+1] + da*dx[i+1];
        dy[i+2] = dy[i+2] + da*dx[i+2];
        dy[i+3] = dy[i+3] + da*dx[i+3];
        }
    }


/*
 Forma o produto escalar de dois vetores.
 Jack Dongarra, linpack, 3/11/78.
 Versão ROLLED
*/
static TIPO_REAL ddot_ur(int n,TIPO_REAL *dx,int incx,TIPO_REAL *dy,int incy)

    {
    TIPO_REAL dtemp;
    int i,ix,iy,m;

    dtemp = ZERO;

    if (n <= 0)
        return(ZERO);

    if (incx != 1 || incy != 1)
        {

        /* código para incrementos desiguais ou incrementos iguais != 1 */

        ix = 0;
        iy = 0;
        if (incx < 0) ix = (-n+1)*incx;
        if (incy < 0) iy = (-n+1)*incy;
        for (i = 0;i < n; i++)
            {
            dtemp = dtemp + dx[ix]*dy[iy];
            ix = ix + incx;
            iy = iy + incy;
            }
        return(dtemp);
        }

    /* código para ambos os incrementos iguais a 1 */

    m = n % 5;
    if (m != 0)
        {
        for (i = 0; i < m; i++)
            dtemp = dtemp + dx[i]*dy[i];
        if (n < 5)
            return(dtemp);
        }
    for (i = m; i < n; i = i + 5)
        {
        dtemp = dtemp + dx[i]*dy[i] +
        dx[i+1]*dy[i+1] + dx[i+2]*dy[i+2] +
        dx[i+3]*dy[i+3] + dx[i+4]*dy[i+4];
        }
    return(dtemp);
    }


/*
** Escala um vetor por uma constante.
** Jack Dongarra, linpack, 3/11/78.
** UNVersão ROLLED
*/
static void dscal_ur(int n,TIPO_REAL da,TIPO_REAL *dx,int incx)

    {
    int i,m,nincx;

    if (n <= 0)
        return;
    if (incx != 1)
        {

        /* code for increment not equal to 1 */

        nincx = n*incx;
        for (i = 0; i < nincx; i = i + incx)
            dx[i] = da*dx[i];
        return;
        }

    /* code for increment equal to 1 */

    m = n % 5;
    if (m != 0)
        {
        for (i = 0; i < m; i++)
            dx[i] = da*dx[i];
        if (n < 5)
            return;
        }
    for (i = m; i < n; i = i + 5)
        {
        dx[i] = da*dx[i];
        dx[i+1] = da*dx[i+1];
        dx[i+2] = da*dx[i+2];
        dx[i+3] = da*dx[i+3];
        dx[i+4] = da*dx[i+4];
        }
    }


/*
** Encontra o índice do elemento com valor absoluto máximo.
** Jack Dongarra, linpack, 3/11/78.
*/
static int idamax(int n,TIPO_REAL *dx,int incx)

    {
    TIPO_REAL dmax;
    int i, ix, itemp;

    if (n < 1)
        return(-1);
    if (n ==1 )
        return(0);
    if(incx != 1)
        {

        /* código para incremento não igual a 1 */

        ix = 1;
        dmax = fabs((double)dx[0]);
        ix = ix + incx;
        for (i = 1; i < n; i++)
            {
            if(fabs((double)dx[ix]) > dmax)
                {
                itemp = i;
                dmax = fabs((double)dx[ix]);
                }
            ix = ix + incx;
            }
        }
    else
        {

        /* código para incremento igual a 1 */

        itemp = 0;
        dmax = fabs((double)dx[0]);
        for (i = 1; i < n; i++)
            if(fabs((double)dx[i]) > dmax)
                {
                itemp = i;
                dmax = fabs((double)dx[i]);
                }
        }
    return (itemp);
    }


static TIPO_REAL second(void)
{
	return (
				(TIPO_REAL)
				(
					(TIPO_REAL)clock()
					/
					(TIPO_REAL)CLOCKS_PER_SEC
				)
			);
}


