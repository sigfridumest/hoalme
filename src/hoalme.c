#include <math.h>
#include <R.h>

/*
void trace(double *x, double *y, double *w, int *n){
	int i, j, N;
	N = *n;
	double **X, **Y, *z, T = 0;

	X = (double **) R_alloc(N, sizeof(double*));
	Y = (double **) R_alloc(N, sizeof(double*)); //three dimensional array for derivatives
	z = (double *) R_alloc(N, sizeof(double));

//extracting and putting into matrices: comes by rows, one matrix after another one (derivatives)
	j = 0;
	for( i = 0; i < N; i++){
		z[ i ] = 0;
		X[ i ] = x + j;
		Y[ i ] = y + j;
		j = j + N;
	};

	for(i = 0; i < N; i++){
		for(j = 0; j < N; j++){
			z[ i ] += X[ i ][ j ] * Y[ j ][ i ]; 
		};
	};

	for(i = 0; i < N; i++){
		*w += z[ i ];
	};
}
*/

void cholderiv(double *x, double *y, int *n, int *d){

	int i, j, k, l, m, N, D, T;
	N = *n;
	D = *d;
	T = N*(N +1)*.5;

	double **L, ***Lder;

	L = (double **) R_alloc(N, sizeof(double*));
	Lder = (double ***) R_alloc(D, sizeof(double**)); //three dimensional array for derivatives

	for( i = 0; i < D; i++){
		Lder[ i ] = (double **) R_alloc(N, sizeof(double*));
	};


//extracting and putting into matrices: comes by rows, one matrix after another one (derivatives)

	j = 0;
	for( i = 0; i < N; i++){
		L[ i ] = x + j;
		j = j + i + 1;
	};

	for( k = 0; k < D; k++){
		m = k*T;
		j = 0;
		for( i = 0; i < N; i++){
			Lder[ k ][ i ] = y + j + m;
			j = j + i + 1;
		};
	};


//decomposing and deriving

	for(m = 0; m < N; m++){
		//Pivot
		L[ m ][ m ] = sqrt(L[ m ][ m ]);
		for(k = 0; k < D; k++){
			Lder[ k ][ m ][ m ] = .5*Lder[ k ][ m ][ m ]/L[ m ][ m ];
		};
		if( m < (N-1)){
		//adjusting lead column
			for(j = (m+1); j < N; j++)	{
				L[ j ][ m ] = L[ j ][ m ]/L[ m ][ m ];
				for(k = 0; k < D; k++)
					Lder[ k ][ j ][ m ] = (  Lder[ k ][ j ][ m ] - L[ j ][ m ]*Lder[ k ][ m ][ m ]  )/L[ m ][ m ];
			};
		//row operations
			for(j = (m+1); j < N; j++)	{
				for(i = j; i < N; i++){
					L[ i ][ j ] = L[ i ][ j ]-L[ i ][ m ]*L[ j ][ m ];
					for(k = 0; k < D; k++)
						Lder[ k ][ i ][ j ] = Lder[ k ][ i ][ j ] - Lder[ k ][ i ][ m ]*L[ j ][ m ] - L[ i ][ m ]*Lder[ k ][ j ][ m ];
				};
			};
		// ends if
		};
	// ends top for
	};
}
