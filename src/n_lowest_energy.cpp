#include "n_lowest_energy.hpp"

float* NLowestE;
int N;
void create_n_lowest_E(int n)
{
	N = n;
	NLowestE = (float *)malloc(sizeof(float) * N);
	for(int i = 0; i < n; i++)
	{
		NLowestE[i] = 9999999;
	}
}

void destroy_n_lowest_E()
{
	free(NLowestE);
}

void shift_following(float E, int i)
{
	float prev = NLowestE[i];
	float next;
	NLowestE[i] = E;
	for(int j = i + 1; j < N - 1; j++)
	{
		next = NLowestE[j];
		NLowestE[j] = prev;
		prev = next;
	}
}

void print_lowest()
{
	printf("Lowest List:\n");
	for(int i = 0; i < N; i++)
	{
		printf("[%d] = %f\n", i, NLowestE[i]);
	}
}

int add_to_n_lowest_E(float E)
{
	int i = 0;
	if(E < NLowestE[0])
	{
		//printf("Replacing [%d] = %f with: %f\n", i, NLowestE[i], E);
		shift_following(E, 0);
		//print_lowest();
		return 0;
	}
	for(i = N-1; i >= 1; i--)
	{
		if(E > NLowestE[i])
		{
			break;	
		}
	}
	i++;
	if(i <= N-1)
	{
		//printf("Replacing [%d] = %f with: %f\n", i, NLowestE[i], E);
		shift_following(E, i);
		//print_lowest();
		return i;
	}
	return -1;
}
