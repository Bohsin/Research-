#define _CRT_RAND_S

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <memory.h>
#include <malloc.h>
#include <math.h>



struct kiss_state {
	/*
	* Seed variables.  Note that kiss requires a moderately complex
	* seeding using a "seed rng" that we will arbitrarily set to be
	* mt19937_1999 from the GSL.
	*/
	unsigned int x;
	unsigned int y;
	unsigned int z;
	unsigned int c;
};


unsigned int kiss_get(struct kiss_state *vstate);
double kiss_get_double(struct kiss_state *vstate);
double kiss_get_double_better(struct kiss_state *vstate);
void init_kiss(struct kiss_state *vstate);
unsigned int randInt();
int kiss_int(struct kiss_state* random_seed, int n);


int main()
{
	int i,d,p,nline,num;
	double b,c;
	int r[10] = {0,0,0,0,0,0,0,0,0,0};
	double sum=0.0,ave=0.0;
	int niteration = 100;
	unsigned int u;


	struct kiss_state* vstate= (struct kiss_state*)malloc(sizeof(struct kiss_state));
	vstate->x = 123456789;
	vstate->y = 987654321;
	vstate->z = 43219876;
	vstate->c = 6543217;

	printf("\n x=%u, y=%u, z=%u, c=%u \n", vstate->x, vstate->y, vstate->z, vstate->c);

	printf("\n UINT_MAX=%lu, ULONG_MAX=%llu, ULLONG_MAX=%llu \n", UINT_MAX, ULONG_MAX, ULLONG_MAX);


	init_kiss(vstate);
	printf("\n %d outputs of kiss_get()%5\n", niteration);
	for (i = 0; i < 500; i++)
	{
		printf("%u ", kiss_get(vstate)%5);
		if (i % 20 == 19) printf("\n");
	}

	init_kiss(vstate);
	printf("\n x=%u, y=%u, z=%u, c=%u \n", vstate->x, vstate->y, vstate->z, vstate->c);

	printf("\n %d outputs of kiss_get_double_better()\n",niteration);
	for (i = 0; i < niteration; i++)
	{
		b = kiss_get_double_better(vstate);
		printf("%f ", b);
		if (i % 20 == 19) printf("\n");
		if (b > 1.00)
		{
			printf("\n error \n");
			break;
		}
		
		d = round(b*1000);
		p = d / 100;
		switch (p)
		{
		case 0:
			r[0] += 1;
			break;
		case 1:
			r[1] += 1;
			break;
		case 2:
			r[2] += 1;
			break;
		case 3:
			r[3] += 1;
			break;
		case 4:
			r[4] += 1;
			break;
		case 5:
			r[5] += 1;
			break;
		case 6:
			r[6] += 1;
			break;
		case 7:
			r[7] += 1;
			break;
		case 8:
			r[8] += 1;
			break;
		case 9:
			r[9] += 1;
			break;
		default:
			break;
		}
		sum += b;
	}

	ave = sum / niteration;

	printf("\n Histogram of random real number between 0 and 1: \n");
	for (i = 0; i < 10; i++)
	{
		if (i<9) 
			printf("Number between 0.%d and 0.%d : %d\n",i,i+1,r[i]);
		else
			printf("Number between 0.%d and 1.0: %d\n", i, r[i]);
	}
	printf("\n Average is: %f\n", ave);

	free(vstate);
	getchar();
	return 0;
}




/*JKISS RNG */
unsigned int kiss_get(struct kiss_state *vstate)
{

	struct kiss_state *state = vstate;
	unsigned long long t;
	state->x = 314527869 * state->x + 1234567;
	state->y ^= state->y << 5;
	state->y ^= state->y >> 7;
	state->y ^= state->y << 22;
	t = 4294584393ULL * state->z + state->c;
	state->c = t >> 32;
	state->z = t;
	unsigned int a = state->x + state->y + state->z;
	return a;

}


/*A simple way to get random real number between 0 and 1;*/
double kiss_get_double(struct kiss_state *vstate)
{
	double a = UINT_MAX+1.0;
	return kiss_get(vstate) / a;
}

/*A better way to get random real number between 0 and 1;*/
double kiss_get_double_better(struct kiss_state *vstate)
{
	double x;
	unsigned int a, b;
	a = kiss_get(vstate) >> 6; /* Upper 26 bits */
	b = kiss_get(vstate) >> 5; /* Upper 27 bits */
	x = (a * 134217728.0 + b) / 9007199254740992.0;
	return x;
}

/*Initialize KISS with randInd()*/
void init_kiss(struct kiss_state *vstate)
{
	struct kiss_state *state = vstate;
	state->x = randInt();
	while (!(state->y = randInt()))
		printf("error"); /* y must not be zero! */
	state->z = randInt();
	/* We don¡¯t really need to set c as well but let's anyway¡­ */
	/* NOTE: offset c by 1 to avoid z=c=0 */
	state->c = randInt() % 698769068 + 1; /* Should be less than 698769069 */
}

/*Get randome number with rand_s(), in order to seed KISS*/
unsigned int randInt()
{
	unsigned int u,a;
	unsigned int min = 1;
	unsigned int max = 4294967295;
	rand_s(&u);
	a = ((double)u / ((double)(UINT_MAX) + 1.0))*(max - min) + min;
	return a;
}

/*To generate a random integer number between 0 and n-1*/
int kiss_int(struct kiss_state* random_seed, int n)
{
	init_kiss(random_seed);
	int random_int = kiss_get(random_seed) % n;
	return random_int;
}
