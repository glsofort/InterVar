def sum_of_list(list):
    sum=0
    for i in list:
        sum=sum+i
    return(sum)


# 1, [0,0,0,0], [0,1,0,0,0,0], [0,0,0,0,0], 0, [0,0,0,0], [0,0,0,0,0,0,0]

def classfyv2(PVS1,PS,PM,PP,BA1,BS,BP):
    BPS=["Pathogenic","Likely pathogenic","Benign","Likely benign","Uncertain significance"]

    BP7 = BP[6]
    BS1 = BS[0]

    PS_sum=sum_of_list(PS)
    PM_sum=sum_of_list(PM)
    PP_sum=sum_of_list(PP)
    BS_sum=sum_of_list(BS)
    BP_sum=sum_of_list(BP)
    
    # Excel logic:
    # IF(
    #     AND(
    #         SUM(C7:C10)>0; SUM(C12:C14)>0
    #     ); 
    #     "VUS";
    #     IF(
    #         SUM(M5:M8;M10)<>0;
    #         N5;
    #         IF(
    #             SUM(F6:F15)<>0;
    #             G5;
    #             IF(
    #                 SUM(I6:I15)<>0;
    #                 J5;
    #                 IF(
    #                     SUM(P6:P8)<>0;
    #                     Q5;
    #                     "VUS"
    #                 )
    #             )
    #         )
    #     )
    # )

    if (PVS1 == 1 or PS_sum > 0 or PM_sum > 0 or PP_sum > 0) and (BA1 == 1 or BS_sum > 0 or BP_sum > 0):
        return (BPS[4]) # Uncertain significance x

    if (
        # M5
        (BP7 == 1 and BP_sum >= 3)
        or
        # M6
        (BA1 == 1)
        or
        # M7
        (BS_sum > 1)
        or
        # M8
        (BS_sum >= 1 and BP_sum >= 3)
        or
        # M10-1
        (BP_sum > 2 and BP7 == 1)
        or 
        # M10-2
        (BS_sum == 1 and BP_sum == 2 and BP7 == 1)
    ):
        return (BPS[2]) # Benign x

    if (
        # F7
        (PVS1 == 1 and PS_sum > 0)
        or
        # F8
        (PVS1 == 1 and PM_sum > 1)
        or 
        # F9
        (PVS1 == 1 and PM_sum > 0 and PP_sum > 0)
        or
        # F10
        (PVS1 == 1 and PP_sum > 1)
        or
        # F11 
        (PS_sum > 1)
        or
        # F13
        (PS_sum > 0 and PM_sum > 2)
        or
        # F14
        (PS_sum > 0 and PM_sum > 1 and PP_sum > 1)
        or
        # F15
        (PS_sum > 0 and PM_sum > 0 and PP_sum > 3)
    ):
        return (BPS[0]) # Pathogenic x

    if (
        # I6
        (PVS1 == 1 and PM_sum == 1)
        or
        # I8
        (PS_sum == 1 and PM_sum < 3 and PM_sum >= 1)
        or
        # I9
        (PS_sum == 1 and PP_sum >= 2)
        or
        # I10
        (PM_sum > 2)
        or
        # I11
        (PM_sum > 1 and PP_sum > 1)
        or
        # I12
        (PM_sum > 0 and PP_sum > 3)
    ):
        return (BPS[1]) # Likely pathogenic x

    if (
        # P6
        (
            (
                (BP7 == 0 or BP_sum < 3) # !(BP7 == 1 and BP_sum >= 3) !M5
                and
                BA1 == 0 # !M6
                and
                BS_sum <= 1 # !M7
            )
            and
            BS1 == 1
            and
            BS_sum >= 1
            and 
            BS_sum < 3
        )
        or
        # P7
        (
            (
                BA1 == 0 # !M6
                and 
                BS_sum <= 1 # !M7
                and
                (BS_sum == 0 or BP_sum < 3) # !M8
                and
                ((BP_sum <= 2 or BP7 == 0) and (BS_sum != 1 or BP_sum != 2 or BP7 == 0)) # !M10
            )
            and
            BS_sum == 1
            and
            BP_sum >= 1
            and
            BP_sum < 3 
        )
        or
        # P8
        (
            (
                BA1 == 0 # !M6
                and 
                BS_sum <= 1 # !M7
                and
                (BS_sum == 0 or BP_sum < 3) # !M8
                and
                ((BP_sum <= 2 or BP7 == 0) and (BS_sum != 1 or BP_sum != 2 or BP7 == 0)) # !M10
            )
            and
            BP_sum >= 2
        )
    ):
        return (BPS[3]) # Likely benign

    return(BPS[4]) # Uncertain significance


cls = classfyv2(1, [0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0], 0, [1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0] )
print(cls)