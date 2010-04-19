
/* compile with -largeArrayDims */
#include <math.h>
#include <mex.h>


void mexFunction( int nout, mxArray *outs[], int nin, const mxArray *ins[] ) {
	int n, m, num_neighs;
	mwIndex *rowindV, *colptrV;
	double *valueV;

	int i, j, k, *maxI;
	mwIndex *rowindE, *colptrE;
	double *valueE, v, maxV, maxE;

	int i1, j1;

	n=mxGetN(ins[0]);
	m=mxGetM(ins[0]);

	colptrV=mxGetJc(ins[0]);
	rowindV=mxGetIr(ins[0]);
	valueV=mxGetPr(ins[0]);

	num_neighs=(int)(*mxGetPr(ins[1]));

	maxI = (int*) mxMalloc( m*sizeof(int) );
	colptrE = (mwIndex*) mxMalloc( (m+1)*sizeof(mwIndex) );
	rowindE = (mwIndex*) mxMalloc( (m*num_neighs)*sizeof(mwIndex) );
	valueE = (double*) mxMalloc( (m*num_neighs)*sizeof(double) );

	for( i=0;i<m;i++ ) {
		colptrE[i]=i*num_neighs;
		maxI[i]=-1;
	}
	colptrE[m]=m*num_neighs;
	/* memset( rowindE,0,(m*num_neighs)*sizeof(size_t) ); */
	for( i=0;i<m*num_neighs;i++ ) {
		valueE[i]=HUGE_VAL;
		rowindE[i]=0;
	}

	for( i=0;i<m;i++ ) {
		for( j=i+1;j<m;j++ ) {

			v=0.0;
			for( k=0;k<n;k++ ) {
				v=v+pow( (valueV[k*m+i]-valueV[k*m+j]), 2.0 );
			}
			v=sqrt(v);

			if( maxI[i] < 0 ) {
				maxV=0.0;
				for( k=0;k<num_neighs;k++ ) {
					 if( maxV < valueE[i*num_neighs+k] ) {
						maxV=valueE[i*num_neighs+k];
						maxI[i]=k;
					}
				}
			} else {
				maxV=valueE[i*num_neighs+maxI[i]];
			}

			if( v < maxV ) {
				valueE[i*num_neighs+maxI[i]]=v;
				rowindE[i*num_neighs+maxI[i]]=j;
				maxI[i]=-1;
			}

			if( maxI[j] < 0 ) {
				maxV=0.0;
				for( k=0;k<num_neighs;k++ ) {
					if( maxV < valueE[j*num_neighs+k] ) {
						maxV=valueE[j*num_neighs+k];
						maxI[j]=k;
					}
				}

			} else {
				maxV=valueE[j*num_neighs+maxI[j]];
			}

			if( v < maxV ) {
				valueE[j*num_neighs+maxI[j]]=v;
				rowindE[j*num_neighs+maxI[j]]=i;
				maxI[j]=-1;
			}
		}
	}

	maxE = 0.0;
	for( i=0; i<m*num_neighs; i++ ) {
		if( maxE == 0 || maxE < valueE[i] ) {
			maxE = valueE[i];
		}
	}

	for( i=0; i<m*num_neighs; i++ ) {
		valueE[i] = exp( - valueE[i]/maxE );
	}

 	outs[0] = mxCreateSparse( m,m,m*num_neighs,mxREAL );
 	mxSetPr( outs[0],valueE );
 	mxSetIr( outs[0],rowindE );
 	mxSetJc( outs[0],colptrE );

	mxFree( maxI );
	return;
}
