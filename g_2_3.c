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

/*It's for flexible array. For this project, it's for unfrozen levels.*/
struct List_int
{
	int length; /*length of the list.*/
	int list[];
};

/*List of unfrozen levels.*/
struct Unfrozen_list
{
	int length1; /*the length of unfrozen levels of chain 1*/
	int length2; /*the length of unfrozen levels of chain 2*/
	int length3; /*the length of unfrozen levels of chain 3*/
	int list[];
};

struct State_Proposal
{
	int energyp;/*energy of proposal state*/
	double P_itop; /*probability of transtion from initial state to proposal state*/
	double P_ptoi;/*probability of transtion from proposal state to initial state*/
	int length1;/*Length of the new unfrozen levels of chain 1*/
	int length2;/*Length of the new unfrozen levels of chain 2*/
	int length3;/*Length of the new unfrozen levels of chain 3*/
	int Microstatep_NewUnfrozenLevel[];/*proposal microstate, with length being 3 times of stateSize calculated in main() and unfrozen levels of proposal state */
};

struct MicroLevel
{
	int chain_num;/*number of chain*/
	int level;
};

struct MicroLevels
{
	int chain_num;
	int length;
	int list[];
};

#define MAXSTEP 1 /*maximum step that a donor can move*/

int VecDotProduct(int *v1, int *v2, int vlength);/*vector dot product, integer type*/
struct Unfrozen_list* UnfrozenLevels(int* Microstate, int stateSize);/*unfrozen levels of state*/
int UnfrozenLevelDeterm_Chain1(int level, int* chain1, int* chain2, int* chain3, int stateSize);/*To determine the 'level' of chain1 is unfrozen or not*/
int UnfrozenLevelDeterm_Chain2(int level, int* chain1, int* chain2, int* chain3, int stateSize);/*To determine the 'level' of chain2 is unfrozen or not*/
int UnfrozenLevelDeterm_Chain3(int level, int* chain1, int* chain2, int* chain3, int stateSize);/*To determine the 'level' of chain3 is unfrozen or not*/
struct List_int* AppendElement(struct List_int* vector, int newElement);/*Add 'newElement' to the end of 'vector'*/
struct List_int* Union_List_int(struct List_int* list1, struct List_int* list2);/*Union of two lists of type struct List_int*/
struct State_Proposal* PDisFun(struct kiss_state* random_seed, int* microstate, int energy, struct Unfrozen_list* unfrozenLevels, int stateSize);
struct MicroLevel* Level_Donor(int nd, int* microstate, struct Unfrozen_list* unfrozenLevel);
struct MicroLevels* Acceptors(int* microstate, struct MicroLevel* level_donor, int stateSize);
struct MicroLevels* Acceptors_chain1(int* microstate, struct MicroLevel* level_donor, int stateSize);
struct MicroLevels* Acceptors_chain2(int* microstate, struct MicroLevel* level_donor, int stateSize);
struct MicroLevels* Acceptors_chain3(int* microstate, struct MicroLevel* level_donor, int stateSize);
struct MicroLevels* AddElement_begin(struct MicroLevels* vector, int newElement);
struct MicroLevels* AddElement_end(struct MicroLevels* vector, int newElement);
struct Unfrozen_list* ReplaceUnfrozenLevels(struct Unfrozen_list* unfrozenLevels1_p, struct Unfrozen_list* unfrozenLevels1);
int* ReplaceMicrostate(int* microState1_p, int* microState1, int stateSize);
double Minimum(double x, double y);

/*Functions of random number generators*/
unsigned int kiss_get(struct kiss_state *vstate);
double kiss_get_double_better(struct kiss_state *vstate);
void init_kiss(struct kiss_state *vstate);
unsigned int randInt();
int kiss_int(struct kiss_state* random_seed, int n);



int main()
{
	int i, j;
	unsigned long k, n;

	double tem = 0.1; /*temperature*/
	int lengthH = (int)(ceil(tem * 5)*2);
	int	grdStLength = lengthH + (int)(ceil(tem)*2);/*length of "Fermi sea of ground state"*/
	int numParTotal = (grdStLength/2) * 3;/*the total number of particles*/
	int EfH = grdStLength - 1;/*Fermi energy*/
	int stateSize = grdStLength + lengthH;/*size of of state*/
	int *energyLevels = (int *)malloc(sizeof(int) *stateSize);
	for (i = 0; i < stateSize; i++) energyLevels[i] = i;/*energy levels, with the energy of lowest level being 0*/

	/*Initial state 1--ground state*/
	int *initlState1 = (int *)malloc(sizeof(int)*stateSize);
	memset(initlState1, 0, sizeof(int)*stateSize);
	/*for (i = 0; i < (grdStLength/2); i++)
	{
		initlState1[2*i] = 2;
		initlState1[2*i+1] = 1;
	}*/
	initlState1[0] = 2;
	initlState1[1] = 1;
	initlState1[2] = 2;
	initlState1[3] = 1;
	initlState1[4] = 2;
	initlState1[5] = 1;
	initlState1[6] = 1;
	initlState1[7] = 1;
	initlState1[8] = 1;
	initlState1[10] = 2;
	initlState1[12] = 1;
	initlState1[15] = 1;
	initlState1[17] = 1;
	initlState1[18] = 1;
	
	int *chain1 = (int *)malloc(sizeof(int)*stateSize);
	memset(chain1, 0, sizeof(int)*stateSize);
	/*for (i = 0; i < (grdStLength / 2); i++)
	{
		chain1[2 * i] = 1;
		chain1[2 * i + 1] = 0;
	}*/
	chain1[0] = 1;
	chain1[2] = 1;
	chain1[4] = 1;
	chain1[6] = 1;
	chain1[10] = 1;
	chain1[15] = 1;
	int *chain2 = (int *)malloc(sizeof(int)*stateSize);
	memset(chain2, 0, sizeof(int)*stateSize);
	/*for (i = 0; i < (grdStLength / 2); i++)
	{
		chain2[2 * i] = 0;
		chain2[2 * i + 1] = 1;
	}*/
	chain2[1] = 1;
	chain2[3] = 1;
	chain2[5] = 1;
	chain2[8] = 1;
	chain2[12] = 1;
	chain2[18] = 1;
	int *chain3 = (int *)malloc(sizeof(int)*stateSize);
	memset(chain3, 0, sizeof(int)*stateSize);
	/*for (i = 0; i < (grdStLength / 2); i++)
	{
		chain3[2 * i] = 1;
		chain3[2 * i + 1] = 0;
	}*/
	chain3[0] = 1;
	chain3[2] = 1;
	chain3[4] = 1;
	chain3[7] = 1;
	chain3[10] = 1;
	chain3[17] = 1;

	int *microState1 = (int *)malloc(sizeof(int)*stateSize * 3);
	memset(microState1, 0, sizeof(int)*stateSize * 3);
	for (i = 0; i < stateSize; i++)
	{
		microState1[i] = chain1[i];
		microState1[stateSize + i] = chain2[i];
		microState1[2 * stateSize + i] = chain3[i];
	}

	int energyInitlState1 = VecDotProduct(initlState1, energyLevels, stateSize);/*Energy of the initial state 1*/
	struct Unfrozen_list* unfrozenLevels1 = UnfrozenLevels(microState1, stateSize);/*unfrozen levels of initial state 1*/

	printf("Temperature is %f\n", tem);
	printf("Size of state: %d\n", stateSize);
	printf("Energy levels:\n");
	for (i = 0; i < stateSize; i++)	printf("%d ", energyLevels[i]);
	printf("\n");
	printf("Initial state 1:\n");
	for (i = 0; i < stateSize; i++) printf("%d ", initlState1[i]);
	printf("\n");
	printf("Microstate:\n");
	printf("Chain 1:\n");
	for (i = 0; i < stateSize; i++) printf("%d ", chain1[i]);
	printf("\n");
	printf("Chain 2:\n");
	for (i = 0; i < stateSize; i++) printf("%d ", chain2[i]);
	printf("\n");
	printf("Chain 3:\n");
	for (i = 0; i < stateSize; i++) printf("%d ", chain3[i]);
	printf("\n");
	printf("Energy of the initial state 1: %d\n", energyInitlState1);
	printf("Unfrozen levels of initial state 1:\n");
	printf("Chain 1:\n");
	for (i = 0; i < (unfrozenLevels1->length1); i++) printf("%d ", unfrozenLevels1->list[i]);
	printf("\n");
	printf("Chain 2:\n");
	for (i = 0; i < (unfrozenLevels1->length2); i++) printf("%d ", unfrozenLevels1->list[unfrozenLevels1->length1 + i]);
	printf("\n");
	printf("Chain 3:\n");
	for (i = 0; i < (unfrozenLevels1->length3); i++) printf("%d ", unfrozenLevels1->list[unfrozenLevels1->length1 + unfrozenLevels1->length2 + i]);
	printf("\n");


	struct kiss_state* random_seed = (struct kiss_state*)malloc(sizeof(struct kiss_state));
	random_seed->x = 123456789;
	random_seed->y = 987654321;
	random_seed->z = 43219876;
	random_seed->c = 6543217;

	unsigned long nIteration[1] = { 1e5 };
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

		for (k = 1; k <= (n - 1); k++)
		{
			rand_0to1 = random_number[k - 1];

			/*printf("Random number is %f\n", rand_0to1);*/

			struct State_Proposal* STATE1_P = PDisFun(random_seed, microState1, energy1, unfrozenLevels1, stateSize);
			int energy1_p = STATE1_P->energyp;
			double Prob1_itop = STATE1_P->P_itop;
			double Prob1_ptoi = STATE1_P->P_ptoi;
			int *microState1_p = (int *)malloc(sizeof(int) * stateSize * 3);
			for (i = 0; i < 3 * stateSize; i++)
				microState1_p[i] = STATE1_P->Microstatep_NewUnfrozenLevel[i];
			int l1 = STATE1_P->length1;
			int l2 = STATE1_P->length2;
			int l3 = STATE1_P->length3;
			struct Unfrozen_list* unfrozenLevels1_p = (struct Unfrozen_list*)malloc(sizeof(struct Unfrozen_list) + (l1 + l2 + l3) * sizeof(int));
			unfrozenLevels1_p->length1 = l1;
			unfrozenLevels1_p->length2 = l2;
			unfrozenLevels1_p->length3 = l3;
			for (i = 0; i < (l1 + l2 + l3); i++)
				unfrozenLevels1_p->list[i] = STATE1_P->Microstatep_NewUnfrozenLevel[3 * stateSize + i];

			double Rmh = exp((energy1 - energy1_p) / tem)*(Prob1_ptoi / Prob1_itop);

			/*printf("\nMetropolis-Hastings ratio is %f\n", Rmh);*/

			if (rand_0to1 < Minimum(1.0, Rmh))
			{
				for (i = 0; i < stateSize; i++)
					stateMatrix1[k][i] = microState1_p[i] + microState1_p[stateSize + i] + microState1_p[2 * stateSize + i];
				energy1 = energy1_p;
				microState1 = ReplaceMicrostate(microState1_p, microState1, stateSize);
				unfrozenLevels1 = ReplaceUnfrozenLevels(unfrozenLevels1_p, unfrozenLevels1);
			}
			else
			{
				for (i = 0; i < stateSize; i++)
					stateMatrix1[k][i] = stateMatrix1[k - 1][i];
			}


			/*printf("Next state is\n");
			for (i = 0; i < stateSize; i++)
			printf("%d ", stateMatrix1[k][i]);
			printf("\n");
			printf("Microstate:\n");
			printf("Chain 1:\n");
			for (i = 0; i < stateSize; i++) printf("%d ", microState1[i]);
			printf("\n");
			printf("Chain 2:\n");
			for (i = 0; i < stateSize; i++) printf("%d ", microState1[stateSize + i]);
			printf("\n");
			printf("Chain 3:\n");
			for (i = 0; i < stateSize; i++) printf("%d ", microState1[2 * stateSize + i]);
			printf("\n");*/


			/*printf("Unfrozen levels of next state are\n");
			printf("Chain 1:\n");
			for (i = 0; i < (unfrozenLevels1->length1); i++) printf("%d ", unfrozenLevels1->list[i]);
			printf("\n");
			printf("Chain 2:\n");
			for (i = 0; i < (unfrozenLevels1->length2); i++) printf("%d ", unfrozenLevels1->list[unfrozenLevels1->length1 + i]);
			printf("\n");
			printf("Chain 3:\n");
			for (i = 0; i < (unfrozenLevels1->length3); i++) printf("%d ", unfrozenLevels1->list[unfrozenLevels1->length1 + unfrozenLevels1->length2 + i]);
			printf("\n");*/

			/*printf("Energy of next state is %d\n", energy1);*/


			free(STATE1_P);
			free(microState1_p);
			free(unfrozenLevels1_p);

		}


		FILE *fp;
		fopen_s(&fp,"data_0.1.txt", "a");
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
	free(chain1);
	free(chain2);
	free(chain3);
	free(microState1);
	free(unfrozenLevels1);
	free(random_seed);

	printf("\nEnd\n");
	getchar();
	return 0;
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


/*To find all the unfrozen levels of 'state'*/
struct Unfrozen_list* UnfrozenLevels(int* Microstate, int stateSize)
{
	int i;
	int *chain1 = (int *)malloc(sizeof(int)*stateSize);
	int *chain2 = (int *)malloc(sizeof(int)*stateSize);
	int *chain3 = (int *)malloc(sizeof(int)*stateSize);
	for (i = 0; i < stateSize; i++)
	{
		chain1[i] = Microstate[i];
		chain2[i] = Microstate[stateSize + i];
		chain3[i] = Microstate[stateSize * 2 + i];
	}
	struct List_int* unfrozenChain1 = (struct List_int*)malloc(sizeof(struct List_int));
	unfrozenChain1->length = 0;
	struct List_int* unfrozenChain2 = (struct List_int*)malloc(sizeof(struct List_int));
	unfrozenChain2->length = 0;
	struct List_int* unfrozenChain3 = (struct List_int*)malloc(sizeof(struct List_int));
	unfrozenChain3->length = 0;

	for (i = 2; i < stateSize; i++)
	{
		if (UnfrozenLevelDeterm_Chain1(i, chain1, chain2, chain3, stateSize) == 1)
			unfrozenChain1 = AppendElement(unfrozenChain1, i);
		if (UnfrozenLevelDeterm_Chain2(i, chain1, chain2, chain3, stateSize) == 1)
			unfrozenChain2 = AppendElement(unfrozenChain2, i);
		if (UnfrozenLevelDeterm_Chain3(i, chain1, chain2, chain3, stateSize) == 1)
			unfrozenChain3 = AppendElement(unfrozenChain3, i);
	}

	struct List_int* union_unfrozenchain12 = Union_List_int(unfrozenChain1, unfrozenChain2);
	struct List_int* union_unfrozenchain123 = Union_List_int(union_unfrozenchain12, unfrozenChain3);

	int l1 = unfrozenChain1->length;
	int l2 = unfrozenChain2->length;
	int l3 = unfrozenChain3->length;
	int l = l1 + l2 + l3;
	struct Unfrozen_list* unfrozenList = (struct Unfrozen_list*)malloc(sizeof(struct Unfrozen_list) + sizeof(int)*l);
	unfrozenList->length1 = l1;
	unfrozenList->length2 = l2;
	unfrozenList->length3 = l3;
	for (i = 0; i < l; i++)
		unfrozenList->list[i] = union_unfrozenchain123->list[i];
	
	free(chain1);
	free(chain2);
	free(chain3);
	free(unfrozenChain1);
	free(unfrozenChain2);
	free(unfrozenChain3);
	free(union_unfrozenchain12);
	free(union_unfrozenchain123);

	return unfrozenList;
}


/*To determine the 'level' of chain1 is unfrozen or not*/
int UnfrozenLevelDeterm_Chain1(int level, int* chain1, int* chain2, int* chain3, int stateSize)
{
	int flag = 0;
	if (chain1[level]==1)
	{
		if (chain2[level-1]==0)
			flag = 1;
		if (flag == 0 && level < (stateSize - 1))
		{
			if (chain3[level] == 0)
				flag = 1;
		}
	}
	return flag;
}

/*To determine the 'level' of chain2 is unfrozen or not*/
int UnfrozenLevelDeterm_Chain2(int level, int* chain1, int* chain2, int* chain3, int stateSize)
{
	int flag = 0;
	if (chain2[level] == 1)
	{
		if (chain3[level - 1] == 0)
			flag = 1;
		if (flag == 0 && level < (stateSize - 1))
		{
			if (chain1[level+1] == 0)
				flag = 1;
		}
	}
	return flag;
}

/*To determine the 'level' of chain3 is unfrozen or not*/
int UnfrozenLevelDeterm_Chain3(int level, int* chain1, int* chain2, int* chain3, int stateSize)
{
	int flag = 0;
	if (chain3[level] == 1)
	{
		if (chain1[level] == 0)
			flag = 1;
		if (flag == 0 && level < (stateSize - 1))
		{
			if (chain2[level + 1] == 0)
				flag = 1;
		}
	}
	return flag;
}

/*Add 'newElement' to the end of 'vector'*/
struct List_int* AppendElement(struct List_int* vector, int newElement)
{

	int i;
	int len = vector->length;
	struct List_int* newVector = (struct List_int*)malloc(sizeof(struct List_int) + (len + 1) * sizeof(int));
	newVector->length = len + 1;
	for (i = 0; i < len; i++)
		newVector->list[i] = vector->list[i];
	newVector->list[len] = newElement;
	free(vector);
	return newVector;
}



/*Union of two lists of type struct List_int*/
struct List_int* Union_List_int(struct List_int* list1, struct List_int* list2)
{
	int i;
	int l1 = list1->length;
	int l2 = list2->length;
	int l = l1 + l2;
	struct List_int* union_listint = (struct List_int*)malloc(sizeof(struct List_int) + sizeof(int)*l);
	union_listint->length = l;
	for (i = 0; i < l1; i++)
		union_listint->list[i] = list1->list[i];
	for (i = 0; i < l2; i++)
		union_listint->list[l1+i] = list2->list[i];
	return union_listint;
}


/*Proposal distribution function*/
struct State_Proposal* PDisFun(struct kiss_state* random_seed, int* microstate, int energy, struct Unfrozen_list* unfrozenLevels, int stateSize)
{
	int i;

	/*to get level of donor and acceptor*/
	int size_donor = (unfrozenLevels->length1) + (unfrozenLevels->length2) + (unfrozenLevels->length3);
	int size_acceptor = 0;
	int nd = kiss_int(random_seed, size_donor);
	struct MicroLevel* levelDonor = Level_Donor(nd, microstate, unfrozenLevels);
	int level_donor = levelDonor->level;
	int chainNumDA = levelDonor->chain_num;

	/*printf("Donor level at chain %d is %d\n", chainNumDA, level_donor);*/

	struct MicroLevels* accepLevels = Acceptors(microstate, levelDonor, stateSize);
	size_acceptor = accepLevels->length;

	/*printf("Acceptor levels at chain %d are\n", chainNumDA);
	for (i = 0; i < size_acceptor; i++)
	printf("%d ", accepLevels->list[i]);
	printf("\n");*/

	int na = kiss_int(random_seed, size_acceptor);
	int level_acceptor = accepLevels->list[na];
	struct MicroLevel* levelAcceptor = (struct MicroLevel*)malloc(sizeof(struct MicroLevel));
	levelAcceptor->chain_num = chainNumDA;
	levelAcceptor->level = level_acceptor;

	/*printf("Acceptor level at chain %d is %d\n", chainNumDA, level_acceptor);*/

	double Pif = 1.0 / (size_donor*size_acceptor);

	/*printf("Probalibity of going from current state to proposal state is %f\n", Pif);*/

	int energy_p = energy - level_donor + level_acceptor;

	/*proposal state*/
	int* microstate_p = (int*)malloc(3*stateSize * sizeof(int));
	for (i = 0; i < (3 * stateSize); i++)
		microstate_p[i] = microstate[i];
	microstate_p[(chainNumDA - 1)*stateSize + level_donor] = microstate_p[(chainNumDA - 1)*stateSize + level_donor] - 1;
	microstate_p[(chainNumDA - 1)*stateSize + level_acceptor] = microstate_p[(chainNumDA - 1)*stateSize + level_acceptor] + 1;

	/*printf("Proposal microstate:\n");
	printf("chain 1:\n");
	for (i = 0; i < stateSize; i++)
	printf("%d ", microstate_p[i]);
	printf("\n");
	printf("chain 2:\n");
	for (i = 0; i < stateSize; i++)
		printf("%d ", microstate_p[stateSize+i]);
	printf("\n");
	printf("chain 3:\n");
	for (i = 0; i < stateSize; i++)
		printf("%d ", microstate_p[stateSize*2+i]);
	printf("\n");*/

	struct Unfrozen_list* unfrozenLevels_proposal= UnfrozenLevels(microstate_p, stateSize);
	struct MicroLevels* accepLevels_proposal = Acceptors(microstate_p, levelAcceptor, stateSize);
	int size_donor_p = (unfrozenLevels_proposal->length1) + (unfrozenLevels_proposal->length2) + (unfrozenLevels_proposal->length3);
	int size_acceptor_p = accepLevels_proposal->length;
	double Pfi = 1.0 / (size_donor_p*size_acceptor_p);

	/*printf("Unfrozen levels of proposal state:\n");
	printf("Chain 1:\n");
	for (i = 0; i < (unfrozenLevels_proposal->length1); i++) printf("%d ", unfrozenLevels_proposal->list[i]);
	printf("\n");
	printf("Chain 2:\n");
	for (i = 0; i < (unfrozenLevels_proposal->length2); i++) printf("%d ", unfrozenLevels_proposal->list[unfrozenLevels_proposal->length1 + i]);
	printf("\n");
	printf("Chain 3:\n");
	for (i = 0; i < (unfrozenLevels_proposal->length3); i++) printf("%d ", unfrozenLevels_proposal->list[unfrozenLevels_proposal->length1 + unfrozenLevels_proposal->length2 + i]);
	printf("\n");*/

	/*printf("Acceptor levels of proposal state at chain %d are\n", chainNumDA);
	for (i = 0; i < size_acceptor_p; i++)
	printf("%d ", accepLevels_proposal->list[i]);
	printf("\n");*/

	/*printf("\nProbalibity of going from proposal state to current state is %f\n", Pfi);*/

	int n = 3*stateSize + size_donor_p;
	struct State_Proposal* STATE_P = (struct State_Proposal*)malloc(sizeof(struct State_Proposal) + n * sizeof(int));
	STATE_P->energyp = energy_p;
	STATE_P->P_itop = Pif;
	STATE_P->P_ptoi = Pfi;
	STATE_P->length1 = unfrozenLevels_proposal->length1;
	STATE_P->length2 = unfrozenLevels_proposal->length2;
	STATE_P->length3 = unfrozenLevels_proposal->length3;
	for (i = 0; i < (3 * stateSize); i++)
		STATE_P->Microstatep_NewUnfrozenLevel[i] = microstate_p[i];
	for (i = 0; i < size_donor_p; i++)
		STATE_P->Microstatep_NewUnfrozenLevel[3 * stateSize + i] = unfrozenLevels_proposal->list[i];

	free(levelDonor);
	free(accepLevels);
	free(levelAcceptor);
	free(microstate_p);
	free(unfrozenLevels_proposal);
	free(accepLevels_proposal);

	return STATE_P;
}


/*To find the level of donor*/
struct MicroLevel* Level_Donor(int nd, int* microstate, struct Unfrozen_list* unfrozenLevel)
{
	int l1 = unfrozenLevel->length1;
	int l2 = unfrozenLevel->length2;
	int l3 = unfrozenLevel->length3;
	struct MicroLevel* level_donor = (struct MicroLevel*)malloc(sizeof(struct MicroLevel));
	level_donor->level = unfrozenLevel->list[nd];
	if (nd < l1)
		level_donor->chain_num = 1;
	else
		if (nd < (l1 + l2))
			level_donor->chain_num = 2;
		else
			level_donor->chain_num = 3;
	return level_donor;
}

/*To find all the levels of acceptor*/
struct MicroLevels* Acceptors(int* microstate, struct MicroLevel* level_donor, int stateSize)
{
	int chainNum = level_donor->chain_num;
	struct MicroLevels* accepLevels=NULL;
	switch (chainNum)
	{
	case 1:
		accepLevels = Acceptors_chain1(microstate, level_donor, stateSize);
		break;
	case 2:
		accepLevels = Acceptors_chain2(microstate, level_donor, stateSize);
		break;
	case 3:
		accepLevels = Acceptors_chain3(microstate, level_donor, stateSize);
		break;
	default:
		printf("\nerror\n");
		break;
	}
	return accepLevels;
}

/*To find levels of acceptor for chain 1*/
struct MicroLevels* Acceptors_chain1(int* microstate, struct MicroLevel* level_donor, int stateSize)
{
	int i;
	int *chain1 = (int *)malloc(sizeof(int)*stateSize);
	int *chain2 = (int *)malloc(sizeof(int)*stateSize);
	int *chain3 = (int *)malloc(sizeof(int)*stateSize);
	for (i = 0; i < stateSize; i++)
	{
		chain1[i] = microstate[i];
		chain2[i] = microstate[stateSize + i];
		chain3[i] = microstate[stateSize * 2 + i];
	}
	int flag_down = 1;
	int flag_up = 1;
	int leveld = level_donor->level;
	struct MicroLevels* accepLevels = (struct MicroLevels*)malloc(sizeof(struct MicroLevels));
	accepLevels->chain_num = 1;
	accepLevels->length = 0;
	for (i = 1; i <=MAXSTEP; i++)
	{
		if ((leveld - i) >= 0 && flag_down == 1)
		{
			if (chain2[leveld - i] == 0)
				accepLevels = AddElement_begin(accepLevels, leveld - i);
			else
				flag_down = 0;
		}
		if ((leveld + i) < stateSize && flag_up == 1)
		{
			if (chain3[leveld + i-1] == 0)
				accepLevels = AddElement_end(accepLevels, leveld + i);
			else
				flag_up = 0;
		}
	}
	free(chain1);
	free(chain2);
	free(chain3);
	return accepLevels;
}

/*To find levels of acceptor for chain 2*/
struct MicroLevels* Acceptors_chain2(int* microstate, struct MicroLevel* level_donor, int stateSize)
{
	int i;
	int *chain1 = (int *)malloc(sizeof(int)*stateSize);
	int *chain2 = (int *)malloc(sizeof(int)*stateSize);
	int *chain3 = (int *)malloc(sizeof(int)*stateSize);
	for (i = 0; i < stateSize; i++)
	{
		chain1[i] = microstate[i];
		chain2[i] = microstate[stateSize + i];
		chain3[i] = microstate[stateSize * 2 + i];
	}
	int flag_down = 1;
	int flag_up = 1;
	int leveld = level_donor->level;
	struct MicroLevels* accepLevels = (struct MicroLevels*)malloc(sizeof(struct MicroLevels));
	accepLevels->chain_num = 2;
	accepLevels->length = 0;
	for (i = 1; i <= MAXSTEP; i++)
	{
		if ((leveld - i) >= 0 && flag_down == 1)
		{
			if (chain3[leveld - i] == 0)
				accepLevels = AddElement_begin(accepLevels, leveld - i);
			else
				flag_down = 0;
		}
		if ((leveld + i) < stateSize && flag_up == 1)
		{
			if (chain1[leveld + i] == 0)
				accepLevels = AddElement_end(accepLevels, leveld + i);
			else
				flag_up = 0;
		}
	}
	free(chain1);
	free(chain2);
	free(chain3);
	return accepLevels;
}

/*To find levels of acceptor for chain 3*/
struct MicroLevels* Acceptors_chain3(int* microstate, struct MicroLevel* level_donor, int stateSize)
{
	int i;
	int *chain1 = (int *)malloc(sizeof(int)*stateSize);
	int *chain2 = (int *)malloc(sizeof(int)*stateSize);
	int *chain3 = (int *)malloc(sizeof(int)*stateSize);
	for (i = 0; i < stateSize; i++)
	{
		chain1[i] = microstate[i];
		chain2[i] = microstate[stateSize + i];
		chain3[i] = microstate[stateSize * 2 + i];
	}
	int flag_down = 1;
	int flag_up = 1;
	int leveld = level_donor->level;
	struct MicroLevels* accepLevels = (struct MicroLevels*)malloc(sizeof(struct MicroLevels));
	accepLevels->chain_num = 3;
	accepLevels->length = 0;
	for (i = 1; i <= MAXSTEP; i++)
	{
		if ((leveld - i) >= 0 && flag_down == 1)
		{
			if (chain1[leveld - i + 1] == 0)
				accepLevels = AddElement_begin(accepLevels, leveld - i);
			else
				flag_down = 0;
		}
		if ((leveld + i) < stateSize && flag_up == 1)
		{
			if (chain2[leveld + i] == 0)
				accepLevels = AddElement_end(accepLevels, leveld + i);
			else
				flag_up = 0;
		}
	}
	free(chain1);
	free(chain2);
	free(chain3);
	return accepLevels;
}

/*Add 'newElement' to the beginning of 'vector' of type struct MicroLevels*/
struct MicroLevels* AddElement_begin(struct MicroLevels* vector, int newElement)
{

	int i;
	int len = vector->length;
	int chainNum = vector->chain_num;
	struct MicroLevels* newVector = (struct MicroLevels*)malloc(sizeof(struct MicroLevels) + (len + 1) * sizeof(int));
	newVector->chain_num = chainNum;
	newVector->length = len + 1;
	newVector->list[0] = newElement;
	for (i = 0; i < len; i++)
		newVector->list[i+1] = vector->list[i];
	free(vector);
	return newVector;
}


/*Add 'newElement' to the end of 'vector' of type struct MicroLevels*/
struct MicroLevels* AddElement_end(struct MicroLevels* vector, int newElement)
{

	int i;
	int len = vector->length;
	int chainNum = vector->chain_num;
	struct MicroLevels* newVector = (struct MicroLevels*)malloc(sizeof(struct MicroLevels) + (len + 1) * sizeof(int));
	newVector->chain_num = chainNum;
	newVector->length = len + 1;
	for (i = 0; i < len; i++)
		newVector->list[i] = vector->list[i];
	newVector->list[len] = newElement;
	free(vector);
	return newVector;
}

/*To replace unfronzen levels with those of proposal state*/
struct Unfrozen_list* ReplaceUnfrozenLevels(struct Unfrozen_list* unfrozenLevels1_p, struct Unfrozen_list* unfrozenLevels1)
{
	int i;
	int l1 = unfrozenLevels1_p->length1;
	int l2 = unfrozenLevels1_p->length2;
	int l3 = unfrozenLevels1_p->length3;
	struct Unfrozen_list* newUnfrozenLevels = (struct Unfrozen_list*)malloc(sizeof(struct Unfrozen_list) + (l1 + l2 + l3) * sizeof(int));
	newUnfrozenLevels->length1 = l1;
	newUnfrozenLevels->length2 = l2;
	newUnfrozenLevels->length3 = l3;
	for (i = 0; i < (l1 + l2 + l3); i++)
		newUnfrozenLevels->list[i] = unfrozenLevels1_p->list[i];
	free(unfrozenLevels1);
	return newUnfrozenLevels;
}


/*To replace the microstate with that of proposal state*/
int* ReplaceMicrostate(int* microState1_p, int* microState1, int stateSize)
{
	int i;
	int* newMicrostate = (int*)malloc(3 * stateSize * sizeof(int));
	for (i = 0; i < (3 * stateSize); i++)
		newMicrostate[i] = microState1_p[i];
	free(microState1);
	return newMicrostate;
}

/*To find the minimum of x and y*/
double Minimum(double x, double y)
{
	if (x <= y)
		return x;
	else
		return y;
}




/*Below are the functions of random number generator.*/
/*JKISS RNG */
unsigned int kiss_get(struct kiss_state *state)
{

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



