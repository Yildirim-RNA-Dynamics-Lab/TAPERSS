/*
if(true)
    {
        for(int i = 0; i < TableRowCount; i++)
        {
            for(int j = 0; j < TableRowCount; j++)
            {
                InteractionTableSum[i] += gsl_matrix_get(InteractionTable, j, i);
                if(InteractionTableSum[i] >=  4)
                {
                    TMP_END = true;
                }
            }
        }

        printf("IDs:%17s", " ");
        for(int is = 0; is < TableRowCount; is++)
        {
            printf("%4d", is);
        }    
        printf("\n");
        printf("Start:%15s", " ");
        for(int is = 0; is < TableRowCount; is++)
        {
            printf("%4d", InteractionTableSum[is]);
        }            
        printf("\n");
        printf("\n");
        
        int MinSum, MinVal, MinIdx, NumInteractions = 0;
        while((MinSum = array_min_idx_for_energy(InteractionTableSum, TableRowCount)) != -1)
        {
            MinVal = TableRowCount;
            MinIdx = -1;
            for(int i = 0; i < TableRowCount; i++)
            {
                if(InteractionTableSum[MinSum] != 0 && MinSum != i && gsl_matrix_get(InteractionTable, MinSum, i) != 0)
                {
                    if(InteractionTableSum[i] > 0)
                    {
                        if(InteractionTableSum[i] < MinVal)
                        {
                            MinVal = InteractionTableSum[i];
                            MinIdx = i;
                        }
                    }
                }
            }
            if(MinIdx != -1)
            {
                
                InteractionTableSum[MinSum] = 0;
                InteractionTableSum[MinIdx] = 0;
                NumInteractions++;
                printf("i:%4d, MinIdx:%4d::", MinSum, MinIdx);
                for(int is = 0; is < TableRowCount; is++)
                {
                   printf("%4d", InteractionTableSum[is]);
                }            
                printf("\n");
            }
            else
            {
                InteractionTableSum[MinSum] = 0;
            }
        }
        printf("Interactions Found = %d\n", NumInteractions);
        gsl_matrix_set_zero(InteractionTable);
        memset(InteractionTableSum, 0, TableRowCount * sizeof(int));
        energy_ -= NumInteractions;
        printf("------------------\n");
    }
    */