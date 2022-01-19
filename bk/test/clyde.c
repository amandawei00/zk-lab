#include <stdio.h>
#include <stdlib.h>
int func(int *x, int n){
	int s = 0;

	/*for(int i=0; i<n; i++){
		printf("%d ", x[i]);
	}
	printf("n = ");
	printf("%d ", n);*/
	for (int i=0; i<n; i++){
		/*printf("i = ");
		printf("%d ", i);
		printf("\n");

		printf("x[i] = ");
		printf("%d ", x[i]);
		printf("\n");*/

		s += x[i];
		//printf("s = ");
		//printf("%d ", x[i]);
		//printf("\n");
	}
	//printf("s = ");
	//printf("%d ", s);
	return s;
}

double func2(int n, double *xx){ //xx=[x,y], x*x + 2*y
	return(xx[0] * xx[0] + 2 * xx[1]);
}

double func3(double x){
	return(x + 8);
}

/*
int main(){
	int x[] = {0, 1, 2, 3, 4};
	int n = 5;

	printf("%d ", func(x, n));
	return(0);
}*/
