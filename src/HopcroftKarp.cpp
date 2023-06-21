#include "HopcroftKarp.hpp"

HKQueue HKQ;
size_t SizeV, SizeU;
int32_t* Pair_V; 
int32_t* Pair_U;
int32_t* Layer;
int32_t* HKMemBlock;

enum HKENUM {HK_NIL = 0, HK_INF = -2};

/**  
 	*	Initializes memory for HK arrays and queue structure. Memory is kept contiguous using HKMemblock.
 	*	size_t N_U: Number of vertices in U side of graph (AKA number of positively charged groups)
 	* size_t N_V: Number of vertices in V side of graph (AKA number of negatively charged groups)
**/
void HK_create(size_t N_U, size_t N_V)
{
    //(((N_U + N_V) * (N_U + N_V))/4) == Maximum Number of Possible Edges
    uint64_t memsize = (((N_V + N_U + 2) * 2) + (((N_U + N_V) * (N_U + N_V))/4)) * sizeof(int32_t);
    size_t offset = 0;
    SizeV = N_V + 1;
    SizeU = N_U + 1;
    HKMemBlock = (int32_t*)malloc(memsize);
    HKQ.Initialize((((N_U + N_V) * (N_U + N_V))/4), HKMemBlock, offset);
    offset += ((((N_U + N_V) * (N_U + N_V))/4));
    Pair_V = &HKMemBlock[offset];
    offset += (SizeV);
    Pair_U = &HKMemBlock[offset];
    offset += (SizeU);
    Layer  = &HKMemBlock[offset];
}

void HK_destroy()
{
    free(HKMemBlock);
}

bool HK_BFS(gsl_matrix* Graph)
{
    for(uint i = 1; i < SizeU; i++)
    {
        if(Pair_U[i] == HKENUM::HK_NIL)
        {
            Layer[i] = 0;
            HKQ.Enqueue(i);
        }
        else
        {
            Layer[i] = HKENUM::HK_INF;
        }
    }

    //HKQ.Print();

    Layer[HKENUM::HK_NIL] = HKENUM::HK_INF;

    while(!HKQ.Empty())
    {
        uint32_t U = HKQ.Dequeue();
        if(Layer[U] != Layer[HKENUM::HK_NIL])
        {
            for(uint i = 0; i < Graph->size2; i++)
            {
                if(gsl_matrix_get(Graph, U - 1, i) == 1)
                {
                    uint32_t V = i + 1;
                    if(Layer[Pair_V[V]] == HKENUM::HK_INF)
                    {
                        Layer[Pair_V[V]] = Layer[U] + 1;
                        HKQ.Enqueue(Pair_V[V]);
                    }
                }
            }
        }
    }

    return (Layer[HKENUM::HK_NIL] != HKENUM::HK_INF);
}

bool HK_DFS(int32_t U, gsl_matrix* Graph)
{
    if(U != HKENUM::HK_NIL)
    {
        for(uint i = 0; i < Graph->size2; i++)
        {
            if(gsl_matrix_get(Graph, U - 1, i) == 1)
            {
                uint32_t V = i + 1;
                if(Layer[Pair_V[V]] == Layer[U] + 1)
                {
                    if(HK_DFS(Pair_V[V], Graph))
                    {
                        Pair_V[V] = U;
                        Pair_U[U] = V;
                        return true;
                    }
                }
            }
        }
        Layer[U] = HKENUM::HK_INF;
        return false;
    }//0-120-135-136-23
    return true;
}

uint32_t HK_GetMaxMatching(gsl_matrix* Graph)
{
    uint32_t Matching = 0;

    HKQ.Reset();

    memset(Pair_V, HKENUM::HK_NIL, sizeof(int32_t) * SizeV);
    memset(Pair_U, HKENUM::HK_NIL, sizeof(int32_t) * SizeU);

    //printf("#########New HK Started###########\n");

    while(HK_BFS(Graph))
    {
        for(uint i = 1; i < SizeU; i++)
        {
            if(Pair_U[i] == HKENUM::HK_NIL && HK_DFS(i, Graph))
            {
                Matching++;
            }
        }
    }
    return Matching;
}
