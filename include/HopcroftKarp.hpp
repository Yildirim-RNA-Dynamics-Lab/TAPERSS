#ifndef HOPCROFTKARP_HPP
#define HOPCROFTKARP_HPP

#include "TAPERSS.hpp"

struct HKQueue
{
    int32_t* Queue;
    int32_t  Iterator;
    int32_t  FIFO_Tracker;
    int32_t  IteratorMax;
    int32_t  NumEnqueued;

    void Initialize(uint32_t MaxSize, int32_t* MemBlock, size_t offset)
    {
        Queue = &MemBlock[offset];
        memset(Queue, -1, MaxSize * sizeof(int32_t));
        Iterator = -1;
        FIFO_Tracker = 0;
        IteratorMax = MaxSize - 1;
        NumEnqueued = 0;
    }
    int32_t Dequeue()
    {
        int32_t rtn;
        if(FIFO_Tracker > IteratorMax)
        {
            FIFO_Tracker = 0;
        }
        if(Queue[FIFO_Tracker] == -1)
        {
            printf("Queue Empty: Cannot Dequeue!\n");
            exit(1);
        }
        rtn = Queue[FIFO_Tracker];
        Queue[FIFO_Tracker] = -1;
        FIFO_Tracker++;
        NumEnqueued--;
        return rtn;
    }
    void Enqueue(uint32_t Value)
    {
        Iterator++;
        if(Iterator > IteratorMax)
        {
            if(Queue[0] != -1)
            {
                printf("Queue overloaded!\n");
                exit(1);
            }
            else
            {
                Iterator = 0;
            }
        }
        Queue[Iterator] = Value;
        NumEnqueued++;
    }
    void Reset()
    {
        Iterator = -1;
        FIFO_Tracker = 0;
        NumEnqueued = 0;
    }
    void Print()
    {
        for(int i = 0; i <= IteratorMax; i++)
        {
            printf("[%4d]", Queue[i]);
            if(FIFO_Tracker == i)
            {
                printf("<-FIFO_Tracker");
            }
            if(Iterator == i)
            {
                printf("<-Iterator");
            }
            printf("\n");
        }
        printf("\n");
    }
    bool Empty()
    {
        return (NumEnqueued == 0);
    }
};

void HK_create(size_t N_U, size_t N_V);
void HK_destroy();
uint32_t HK_GetMaxMatching(gsl_matrix* Graph);

#endif
