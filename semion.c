#define _CRT_RAND_S

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <memory.h>
#include <malloc.h>
#include <math.h>


struct kiss_state 
{
	/*Seed variables*/
	unsigned int x;
	unsigned int y;
	unsigned int z;
	unsigned int c;
};

struct State_Proposal
{
	int energyp;/*energy of proposal state*/
	double P_itop; /*probability of transtion from initial state to proposal state*/
	double P_ptoi;/*probability of transtion from proposal state to initial state*/
	int n_NewUnfrozen;/*Length of the new unfrozen levels*/
	int statep_NewUnfrozenLevel[];/*proposal state, with length being the sum of stateSize calculated in main() and unfrozen levels of proposal state */
};

/*It's for flexible array. For this project, it's for unfrozen levels.*/
struct List_int
{
	int length; /*The 1st member is the length of the array.*/
	int list[];
};

#define MSEMION 2 /*m of semion is 2, Haldane statistical parameter g=1/m*/
#define MAXSTEP 1 /*maximum step that a donor can move*/

int VecDotProduct(int *v1, int *v2, int vlength);/*vector dot product, integer type*/
int Unfrozen2Determ(int level, int *state, int stateSize);/*to determine if the 2 is unfrozened level or not*/
struct List_int* UnfrozenLevels(int *state, int stateSize);/*unfrozen levels of state*/
struct List_int* AppendElement(struct List_int* vector, int newElement);
struct State_Proposal* PDisFunSemion(struct kiss_state* random_seed, int* state, int energy, struct List_int* unfrozenLevel, int stateSize);
struct List_int* Acceptors_Semion(int* state, int level_donor, int stateSize);
struct List_int* Acceptors_donor1(int* state, int level_donor, int stateSize);
struct List_int* Acceptors_donor2(int* state, int level_donor, int stateSize);
unsigned int Even1Determ(int* state, int level_donor);
double Minimum(double x, double y);
struct List_int*  ReplaceUnfrozenLevels(struct List_int* unfrozenLevels1_p, struct List_int* unfrozenLevels1);
unsigned int ZeroBetween(int* state, int start, int end);

/*Functions of random number generators*/
unsigned int kiss_get(struct kiss_state *vstate);
double kiss_get_double_better(struct kiss_state *vstate);
void init_kiss(struct kiss_state *vstate);
unsigned int randInt();
int kiss_int(struct kiss_state* random_seed, int n);


int main()
{
	int i,j;
	unsigned long k,n;


	/*parameter initialization*/
	double tem=1; /*temperature*/
	int lengthH = (int)ceil(tem * 10);
	int	grdStLength = lengthH + (int)ceil(2 * tem);/*length of "Fermi sea of ground state"*/
	int numParTotal = grdStLength * 2;/*the total number of particles*/
	int EfH = grdStLength - 1;/*Fermi energy*/
	int stateSize = grdStLength + lengthH;/*size of of state*/
	int *energyLevels = (int *)malloc(sizeof(int) *stateSize);
	for (i = 0; i < stateSize; i++) energyLevels[i] = i;/*energy levels, with the energy of lowest level being 0*/

	/*Initial state 1--ground state*/
	int *initlState1 = (int *)malloc(sizeof(int)*stateSize);
	memset(initlState1, 0, sizeof(int)*stateSize);
	for (i = 0; i < grdStLength; i++) initlState1[i] = MSEMION;
	/*initlState1[0] = 2;
	initlState1[1] = 2;
	initlState1[2] = 2;
	initlState1[3] = 2;
	initlState1[4] = 2;
	initlState1[5] = 2;
	initlState1[6] = 2;
	initlState1[7] = 0;
	initlState1[8] = 0;
	initlState1[9] = 2;
	initlState1[10] = 2;
	initlState1[11] = 2;
	initlState1[12] = 0;
	initlState1[13] = 1;
	initlState1[14] = 1;
	initlState1[15] = 0;
	initlState1[16] = 2;
	initlState1[17] = 0;
	initlState1[18] = 0;
	initlState1[19] = 0;
	initlState1[20] = 0;
	initlState1[21] = 0;*/


	int energyInitlState1 = VecDotProduct(initlState1, energyLevels, stateSize);/*Energy of the initial state 1*/
	struct List_int* unfrozenLevels1 = UnfrozenLevels(initlState1, stateSize);/*unfrozen levels of initial state 1*/


	printf("Temperature is %f\n", tem);
	printf("Size of state: %d\n", stateSize);
	printf("Energy levels:\n");
	for (i = 0; i < stateSize; i++)	printf("%d ", energyLevels[i]);
	printf("\n");
	printf("Initial state 1:\n");
	for (i = 0; i < stateSize; i++) printf("%d ", initlState1[i]);
	printf("\n");
	printf("Energy of the initial state 1: %d\n", energyInitlState1);
	printf("Unfrozen levels of initial state 1:\n");
	for (i = 0; i < (int)(unfrozenLevels1->length); i++)
		printf("%d ", unfrozenLevels1->list[i]);
	printf("\n");

	
	struct kiss_state* random_seed = (struct kiss_state*)malloc(sizeof(struct kiss_state));
	random_seed->x = 123456789;
	random_seed->y = 987654321;
	random_seed->z = 43219876;
	random_seed->c = 6543217;


	unsigned long nIteration[1] = {1e6};
	for (j = 0; j < 1; j++)
	{
		n = nIteration[j];


		/*Initialize a state matrix with 0s, with size n by stateSize*/
		int **stateMatrix1;
		stateMatrix1 = (int**)malloc(n*(sizeof(int*)));
		for (k = 0; k<n; k++)
			stateMatrix1[k] = (int*)malloc(sizeof(int)*stateSize);
		for (i = 0; i < stateSize; i++)
			stateMatrix1[0][i] = initlState1[i];
		for (k = 1; k < n; k++)
			for (i = 0; i < stateSize; i++)
				stateMatrix1[k][i] = 0;


		/*To produce n random numbers between 0 and 1 with KISS*/
		double* random_number = (double*)malloc(n * sizeof(double));
		init_kiss(random_seed);
		for (k = 0; k < n; k++)
			random_number[k] = kiss_get_double_better(random_seed);

		double rand_0to1;
		int energy1 = energyInitlState1;

		for (k = 1; k <=(n-1); k++)
		{
			rand_0to1 = random_number[k-1];

			/*printf("\nRandom number is %f\n", rand_0to1);*/

			struct State_Proposal* STATE1_P = PDisFunSemion(random_seed, stateMatrix1[k - 1], energy1, unfrozenLevels1, stateSize);
			int energy1_p = STATE1_P->energyp;
			double Prob1_itop = STATE1_P->P_itop;
			double Prob1_ptoi = STATE1_P->P_ptoi;
			int n1_newunfrozen = STATE1_P->n_NewUnfrozen;
			int* state1_p = (int*)malloc(stateSize * sizeof(int));
			for (i = 0; i < stateSize; i++)
				state1_p[i] = STATE1_P->statep_NewUnfrozenLevel[i];
			struct List_int* unfrozenLevels1_p = (struct List_int*)malloc(sizeof(struct List_int)+n1_newunfrozen * sizeof(int));
			unfrozenLevels1_p->length = n1_newunfrozen;
			for (i = 0; i < n1_newunfrozen; i++)
				unfrozenLevels1_p->list[i] = STATE1_P->statep_NewUnfrozenLevel[i+stateSize];

			double Rmh = exp((energy1 - energy1_p) / tem)*(Prob1_ptoi / Prob1_itop);

			/*printf("\nMetropolis-Hastings ratio is %f\n", Rmh);*/

			if ( rand_0to1 < Minimum(1.0, Rmh) )
			{
				for (i = 0; i < stateSize; i++)
					stateMatrix1[k][i] = state1_p[i];
				energy1 = energy1_p;
				unfrozenLevels1 = ReplaceUnfrozenLevels(unfrozenLevels1_p, unfrozenLevels1);
			}
			else
			{
				for (i = 0; i < stateSize; i++)
					stateMatrix1[k][i] = stateMatrix1[k-1][i];
			}

			/*printf("\nNext state is\n");
			for (i = 0; i < stateSize; i++)
				printf("%d ", stateMatrix1[k][i]);
			printf("\n");*/


			/*printf("\nUnfrozen levels of next state are\n");
			for (i = 0; i < unfrozenLevels1->length; i++)
				printf("%d ", unfrozenLevels1->list[i]);
			printf("\n");*/

			/*printf("\nEnergy of next state is %d\n", energy1);*/

			free(state1_p);
			free(STATE1_P);
			free(unfrozenLevels1_p);
		}

		

		FILE *fp;
		fopen_s(&fp,"data.txt", "a");
		for (k = 0; k < n; k++)
		{
			for (i = 0; i < stateSize; i++)
			{
				fprintf(fp, "%d ", stateMatrix1[k][i]);
			}
			fprintf(fp, "\n");
		}
		fclose(fp);
		
		for (k = 0; k < n; k++) 
			free(stateMatrix1[k]);
		free(stateMatrix1);
		free(random_number);
	}



	
	free(energyLevels);
	free(initlState1);
	free(unfrozenLevels1);
	free(random_seed);

	printf("\nEnd\n");
	getchar();
	return 0;
}


/*To find all the unfrozen levels of 'state'*/
struct List_int* UnfrozenLevels(int *state, int stateSize)
{
	int i,start2;
	struct List_int* unfrozenLevels = (struct List_int*)malloc(sizeof(struct List_int)+sizeof(int));

	/*To initialize unfrozenLevels*/
	unfrozenLevels->length = 1;
	for (i = 0; i < stateSize; i++)
	{
		if (state[i]>0)
		{
			if (state[i] == 1)
			{
				unfrozenLevels->list[0] = i;
				start2 = i + 1;
				break;
			}
			else
			{
				if (Unfrozen2Determ(i, state, stateSize) == 1)
				{
					unfrozenLevels->list[0] = i;
					start2 = i + 1;
					break;
				}
			}
		}
	}

	/*To include all the remaining unfrozen levels*/
	for (i = start2; i < stateSize; i++)
	{
		if (state[i]>0)
		{
			if (state[i] == 1)
			{
				unfrozenLevels= AppendElement(unfrozenLevels, i);
			}
			else
			{
				if (Unfrozen2Determ(i, state, stateSize) == 1)
				{
					unfrozenLevels = AppendElement(unfrozenLevels, i);
				}
			}
		}
	}
	return unfrozenLevels;
}



/*vector dot product of two int type vectors, i.e. v1 and v2*/
int VecDotProduct(int *v1, int *v2, int vlength)
{
	int i;
	int product = 0;
	for (i = 0; i < vlength; i++)
	{
		product += v1[i] * v2[i];
	}
	return product;
}


/*To determine if a level of 2-particle is unfrozen or not.*/
int Unfrozen2Determ(int level, int *state, int stateSize)
{
	int flag = 0;

	if (level==0)
	{
		if (state[1]==0)
			flag = 1;
	}
	else
	{
		if (level==stateSize-1)
		{
			if (state[level-1]==0)
				flag = 1;
		}
		else
		{
			if (state[level - 1] == 0 || state[level + 1] == 0)
				flag = 1;
		}
	}

	return flag;
}


/*Add 'newElement' to the end 'vector'*/
struct List_int* AppendElement(struct List_int* vector, int newElement)
{
	
	int i;
	int len = vector->length;
	struct List_int* newVector = (struct List_int*)malloc(sizeof(struct List_int) + (len+1)*sizeof(int));
	newVector->length = len + 1;
	for (i = 0; i < len; i++)
		newVector->list[i] = vector->list[i];
	newVector->list[len] = newElement;
	free(vector);
	return newVector;
}





/*Proposal step of MCMC*/
struct State_Proposal* PDisFunSemion(struct kiss_state* random_seed, int* state, int energy, struct List_int* unfrozenLevel, int stateSize)
{
	int i, j;

	/*To get level of donor and acceptor*/
	int size_donor = unfrozenLevel->length;
	int size_acceptor = 0;
	int nd = kiss_int(random_seed, size_donor);
	int level_donor = unfrozenLevel->list[nd];

	/*printf("\nLevel of donor is %d\n", level_donor);*/

	struct List_int* accepLevel = Acceptors_Semion(state, level_donor, stateSize);
	size_acceptor = accepLevel->length;

	/*printf("\nPossible levels of acceptor are\n");
	for (i = 0; i < size_acceptor; i++)
		printf("%d ", accepLevel->list[i]);
	printf("\n");*/

	int na = kiss_int(random_seed, size_acceptor);
	int level_acceptor = accepLevel->list[na];

	/*printf("\nLevel of acceptor is %d\n", level_acceptor);*/

	double Pif = 1.0 / (size_donor*size_acceptor);

	/*printf("\nProbalibity of going from current state to proposal state is %f\n", Pif);*/

	int energy_p = energy - level_donor + level_acceptor;

	/*printf("\nEnergy of proposal state is %d\n", energy_p);*/

	/*proposal state*/
	int* state_p = (int*)malloc(stateSize * sizeof(int));
	
	for (i = 0; i < stateSize; i++)
		state_p[i] = state[i];
	state_p[level_donor] = state_p[level_donor] - 1;
	state_p[level_acceptor] = state_p[level_acceptor] + 1;

	/*printf("\nProposal state is\n");
	for (i = 0; i < stateSize; i++)
		printf("%d ", state_p[i]);
	printf("\n");*/

	struct List_int* unfrozenLevel_proposal = UnfrozenLevels(state_p, stateSize);
    struct List_int* accepLevel_proposal = Acceptors_Semion(state_p, level_acceptor, stateSize);
	int size_donor_p = unfrozenLevel_proposal->length;
	int size_acceptor_p = accepLevel_proposal->length;
	double Pfi = 1.0 / (size_donor_p*size_acceptor_p);

	/*printf("\nUnfrozen levels of proposal state are\n");
	for (i = 0; i < size_donor_p; i++)
		printf("%d ", unfrozenLevel_proposal->list[i]);
	printf("\n");*/

	/*printf("\nAcceptor levels of proposal state are\n");
	for (i = 0; i < size_acceptor_p; i++)
		printf("%d ", accepLevel_proposal->list[i]);
	printf("\n");*/
	
	/*printf("\nProbalibity of going from proposal state to current state is %f\n", Pfi);*/

	int n = stateSize + size_donor_p;
	struct State_Proposal* STATE_P = (struct State_Proposal*)malloc(sizeof(struct State_Proposal) + n * sizeof(int));
	STATE_P->energyp = energy_p;
	STATE_P->P_itop = Pif;
	STATE_P->P_ptoi = Pfi;
	STATE_P->n_NewUnfrozen = size_donor_p;
	for (i = 0; i < stateSize; i++)
		STATE_P->statep_NewUnfrozenLevel[i] = state_p[i];
	j = 0;
	for (i = stateSize; i < n; i++)
	{
		STATE_P->statep_NewUnfrozenLevel[i] = unfrozenLevel_proposal->list[j];
		j++;
	}
		
	free(accepLevel);
	free(state_p);
	free(unfrozenLevel_proposal);
	free(accepLevel_proposal);

	return STATE_P;
}


/*Levels of acceptors around donor 1 or 2*/
struct List_int* Acceptors_Semion(int* state, int level_donor, int stateSize)
{
	int donor = state[level_donor];
	struct List_int* accepLevel;
	if (donor == 2)
		accepLevel = Acceptors_donor2(state, level_donor, stateSize);
	else
		accepLevel = Acceptors_donor1(state, level_donor, stateSize);
	return accepLevel;
}


/*Levels of acceptors around donor 1*/
struct List_int* Acceptors_donor1(int* state, int level_donor, int stateSize)
{
	unsigned int flag_even1 = Even1Determ(state, level_donor);
	int i, start2;/*level*/
	int min = level_donor - MAXSTEP;
	int max = level_donor + MAXSTEP;
	if (min < 0) min = 0;
	if (max > stateSize - 1) max = stateSize - 1;
	int startup = level_donor + 1;
	int startdown = level_donor - 1;
	if (startup > stateSize - 1) startup = stateSize-1;
	if (startdown < 0)startdown = 0;

	struct List_int* accepLevel = (struct List_int*)malloc(sizeof(struct List_int) + sizeof(int));

	/*To initialize accepLevel*/
	unsigned int f_initl = 0;
	accepLevel->length = 1;
	if (flag_even1 == 1)
	{
		/*even 1*/
		for (i =min; i <=startdown; i++)
		{
			if (state[i]==0 && ZeroBetween(state, i, level_donor) == 1)
			{
				accepLevel->list[0] = i;
				start2 = i + 1;
				f_initl = 1;
				break;
			}
			if (state[i]==1 && ZeroBetween(state, i, level_donor)==1)
			{
				accepLevel->list[0] = i;
				start2 = i + 1;
				f_initl = 1;
				break;
			}
		}
		if (f_initl == 0)
		{
			for (i = startup; i <= max; i++)
			{
				if (state[i] == 0 && ZeroBetween(state, level_donor, i)==1)
				{
					accepLevel->list[0] = i;
					start2 = i + 1;
					break;
				}
			}
		}
		
	}
	else
	{
		/*odd 1*/
		for (i = min; i <= startdown; i++)
		{
			if (state[i] == 0 && ZeroBetween(state, i, level_donor) == 1)
			{
				accepLevel->list[0] = i;
				start2 = i + 1;
				f_initl = 1;
				break;
			}

		}
		if (f_initl == 0)
		{
			for (i = startup; i <= max; i++)
			{
				if (state[i] == 0 && ZeroBetween(state, level_donor, i) == 1)
				{
					accepLevel->list[0] = i;
					start2 = i + 1;
					break;
				}
				if (state[i] == 1 && ZeroBetween(state, level_donor, i) == 1)
				{
					accepLevel->list[0] = i;
					start2 = i + 1;
					break;
				}
			}
		}
	}


	/*To include all the remaining accecptor levels*/
	if (flag_even1 == 1)
	{
		/*even 1*/
		if (start2 <= startdown)
		{
			for (i=start2; i<=startdown; i++)
			{
				if (state[i] == 0 && ZeroBetween(state, i, level_donor) == 1)
					accepLevel = AppendElement(accepLevel, i);
				if (state[i] == 1 && ZeroBetween(state, i, level_donor) == 1)
					accepLevel = AppendElement(accepLevel, i);
			}
			for (i = startup; i <= max; i++)
			{
				if (state[i] == 0 && ZeroBetween(state, level_donor, i) == 1)
					accepLevel = AppendElement(accepLevel, i);
			}
		}
		else
		{
			if (start2 == level_donor) start2 = startup;
			if (start2>=startup && start2<=max)
			{
				for (i = start2; i <= max; i++)
				{
					if (state[i] == 0 && ZeroBetween(state, level_donor, i) == 1)
						accepLevel = AppendElement(accepLevel, i);
				}

			}
			
		}
	}
	else
	{
		/*odd 1*/
		if (start2 <= startdown)
		{
			for (i = start2; i <= startdown; i++)
			{
				if (state[i] == 0 && ZeroBetween(state, i, level_donor) == 1)
					accepLevel = AppendElement(accepLevel, i);
			}
			for (i = startup; i <= max; i++)
			{
				if (state[i] == 0 && ZeroBetween(state, level_donor, i) == 1)
					accepLevel = AppendElement(accepLevel, i);
				if (state[i] == 1 && ZeroBetween(state, level_donor, i) == 1)
					accepLevel = AppendElement(accepLevel, i);
			}
		}
		else
		{
			if (start2 == level_donor) start2 = startup;
			if (start2 >= startup && start2 <= max)
			{
				for (i = start2; i <= max; i++)
				{
					if (state[i] == 0 && ZeroBetween(state, level_donor, i) == 1)
						accepLevel = AppendElement(accepLevel, i);
					if (state[i] == 1 && ZeroBetween(state, level_donor, i) == 1)
						accepLevel = AppendElement(accepLevel, i);
				}
			}
		}
	}


	return accepLevel;
}


/*Levels of acceptors around donor 2*/
struct List_int* Acceptors_donor2(int* state, int level_donor, int stateSize)
{
	int i,start2;
	int min = level_donor - MAXSTEP;
	int max = level_donor + MAXSTEP;
	if (min < 0) min = 0;
	if (max > stateSize - 1) max = stateSize - 1;
	struct List_int* accepLevel=(struct List_int*)malloc(sizeof(struct List_int) + sizeof(int));
	
	/*To initialize accepLevel*/
	accepLevel->length = 1;
	for (i = min; i <= max; i++)
	{
		if (i < level_donor)
		{
			if (state[i] == 0 && ZeroBetween(state, i, level_donor) == 1)
			{
				accepLevel->list[0] = i;
				start2 = i + 1;
				break;
			}
		}
		if (i > level_donor)
		{
			if (state[i] == 0 && ZeroBetween(state, level_donor, i) == 1)
			{
				accepLevel->list[0] = i;
				start2 = i + 1;
				break;
			}
		}
		
	}


	/*To include all the remaining accecptor levels*/
	for (i = start2; i <=max; i++)
	{
		if (i < level_donor)
		{
			if (state[i] == 0 && ZeroBetween(state, i, level_donor) == 1)
				accepLevel = AppendElement(accepLevel, i);
		}
		if (i > level_donor)
		{
			if (state[i] == 0 && ZeroBetween(state, level_donor, i) == 1)
				accepLevel = AppendElement(accepLevel, i);
		}
	}

	return accepLevel;
}

/*To determine the 1 is even or odd: return 1 if it's even; return 0 otherwise*/
unsigned int Even1Determ(int* state, int level_donor)
{
	unsigned int flag = 0;
	int sum = 0;
	int i;
	for (i = 0; i <=level_donor; i++)
		sum = sum + state[i];
	if ((sum % 2) == 0) flag = 1;
	return flag;		
}

/*To return 1 if there is no particle between level 'start' and level 'end'; return 0 otherwise*/
unsigned int ZeroBetween(int* state, int start, int end)
{
	int i;
	unsigned int flag = 1;
	if (start<end-1)
	{
		for (i = start+1; i <=end-1; i++)
		{
			if (state[i]!=0)
			{
				flag = 0;
				break;
			}
		}
	}

	return flag;
}


/*To find the minimum of x and y*/
double Minimum(double x, double y)
{
	if (x <= y)
		return x;
	else
		return y;
}



/*To replace unfronzen levels with those of proposal state*/
struct List_int*  ReplaceUnfrozenLevels(struct List_int* unfrozenLevels1_p, struct List_int* unfrozenLevels1)
{
	int i;
	int n = unfrozenLevels1_p->length;
	struct List_int* newUnfrozenLevels = (struct List_int*)malloc(sizeof(struct List_int) + n * sizeof(int));
	newUnfrozenLevels->length = n;
	for (i = 0; i < n; i++)
		newUnfrozenLevels->list[i] = unfrozenLevels1_p->list[i];
	free(unfrozenLevels1);
	return newUnfrozenLevels;
}




/*Below are the functions of random number generator.*/
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
void init_kiss(struct kiss_state* vstate)
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
	unsigned int u, a;
	unsigned int min = 1;
	unsigned int max = 4294967295;
	rand_s(&u);
	a = ((double)u / ((double)(UINT_MAX)+1.0))*(max - min) + min;
	return a;
}

/*To generate a random integer number between 0 and n-1*/
int kiss_int(struct kiss_state* random_seed, int n)
{
	init_kiss(random_seed);
	int random_int = kiss_get(random_seed) % n;
	return random_int;
}

