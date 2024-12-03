def sum_of_list(list):
    sum = 0
    for i in list:
        sum = sum + i
    return sum


# 1, [0,0,0,0], [0,1,0,0,0,0], [0,0,0,0,0], 0, [0,0,0,0], [0,0,0,0,0,0,0]


def classfyv2(PVS1, PS, PM, PP, BA1, BS, BP, PVS):
    BPS = [
        "Pathogenic",
        "Likely pathogenic",
        "Benign",
        "Likely benign",
        "Uncertain significance",
    ]

    BP7 = BP[6]
    BS1 = BS[0]

    PS_sum = sum_of_list(PS)
    PM_sum = sum_of_list(PM)
    PP_sum = sum_of_list(PP)
    BS_sum = sum_of_list(BS)
    BP_sum = sum_of_list(BP)
    PVS_sum = sum_of_list(PVS)

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

    # Exception for only PM2_Supporting and PVS1 exist -> LP
    try:
        if (
            PVS1 == 1
            and PVS_sum == 0
            and PS_sum == 0
            and PM_sum == 0
            and PP_sum == 1
            and BA1 == 0
            and BS_sum == 0
            and BP_sum == 0
        ):
            PM2_Supporting = PP[6]
            if PM2_Supporting == 1:
                return BPS[1]
    except KeyError:
        pass

    # Updating PVS_sum with PVS1 into account
    PVS_sum = PVS_sum + 1 if PVS1 == 1 else PVS_sum

    if (PVS_sum > 0, PS_sum > 0 or PM_sum > 0 or PP_sum > 0) and (
        BA1 == 1 or BS_sum > 0 or BP_sum > 0
    ):
        return BPS[4]  # Uncertain significance x

    # If include PVS1, total PVS class > 2 then it is pathogenic
    if PVS_sum > 1:
        return BPS[0]  # Pathogenic x

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
        return BPS[2]  # Benign x

    if (
        # F7
        (PVS_sum == 1 and PS_sum > 0)
        or
        # F8
        (PVS_sum == 1 and PM_sum > 1)
        or
        # F9
        (PVS_sum == 1 and PM_sum > 0 and PP_sum > 0)
        or
        # F10
        (PVS_sum == 1 and PP_sum > 1)
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
        return BPS[0]  # Pathogenic x

    if (
        # I6
        (PVS_sum == 1 and PM_sum == 1)
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
        return BPS[1]  # Likely pathogenic x

    if (
        # P6
        (
            (
                (BP7 == 0 or BP_sum < 3)  # !(BP7 == 1 and BP_sum >= 3) !M5
                and BA1 == 0  # !M6
                and BS_sum <= 1  # !M7
            )
            and BS1 == 1
            and BS_sum >= 1
            and BS_sum < 3
        )
        or
        # P7
        (
            (
                BA1 == 0  # !M6
                and BS_sum <= 1  # !M7
                and (BS_sum == 0 or BP_sum < 3)  # !M8
                and (
                    (BP_sum <= 2 or BP7 == 0)
                    and (BS_sum != 1 or BP_sum != 2 or BP7 == 0)
                )  # !M10
            )
            and BS_sum == 1
            and BP_sum >= 1
            and BP_sum < 3
        )
        or
        # P8
        (
            (
                BA1 == 0  # !M6
                and BS_sum <= 1  # !M7
                and (BS_sum == 0 or BP_sum < 3)  # !M8
                and (
                    (BP_sum <= 2 or BP7 == 0)
                    and (BS_sum != 1 or BP_sum != 2 or BP7 == 0)
                )  # !M10
            )
            and BP_sum >= 2
        )
    ):
        return BPS[3]  # Likely benign

    return BPS[4]  # Uncertain significance


# Example run
def test():
    evd = {
        # PVS1
        "PVS1": 1,
        # PS
        "PS1": 0,
        "PS2": 0,
        "PS3": 0,
        "PS4": 0,
        "PVS1_Strong": 0,
        "PP3_Strong": 0,
        # PM
        "PM1": 0,
        "PM2": 0,
        "PM3": 0,
        "PM4": 0,
        "PM5": 0,
        "PM6": 0,
        "PVS1_Moderate": 0,
        "PP3_Moderate": 0,
        # PP
        "PP1": 0,
        "PP2": 0,
        "PP3": 0,
        "PP4": 0,
        "PP5": 0,
        "PVS1_Supporting": 0,
        "PM2_Supporting": 0,
        # BA1
        "BA1": 0,
        # BS
        "BS1": 0,
        "BS2": 0,
        "BS3": 0,
        "BS4": 0,
        "BP4_Strong": 0,
        # BP
        "BP1": 0,
        "BP2": 0,
        "BP3": 0,
        "BP4": 0,
        "BP5": 0,
        "BP6": 0,
        "BP7": 0,
        # PVS
        "PM3_VeryStrong": 0,
        "PS2_VeryStrong": 0,
    }

    # List of evidences
    PVS1 = evd["PVS1"]
    PS = [
        evd["PS1"],
        evd["PS2"],
        evd["PS3"],
        evd["PS4"],
        evd["PVS1_Strong"],
        evd["PP3_Strong"],
    ]
    PM = [
        evd["PM1"],
        evd["PM2"],
        evd["PM3"],
        evd["PM4"],
        evd["PM5"],
        evd["PM6"],
        evd["PVS1_Moderate"],
        evd["PP3_Moderate"],
    ]
    PP = [
        evd["PP1"],
        evd["PP2"],
        evd["PP3"],
        evd["PP4"],
        evd["PP5"],
        evd["PVS1_Supporting"],
        evd["PM2_Supporting"],
    ]
    BA1 = evd["BA1"]
    BS = [evd["BS1"], evd["BS2"], evd["BS3"], evd["BS4"], evd["BP4_Strong"]]
    BP = [
        evd["BP1"],
        evd["BP2"],
        evd["BP3"],
        evd["BP4"],
        evd["BP5"],
        evd["BP6"],
        evd["BP7"],
    ]
    PVS = [
        evd["PM3_VeryStrong"],
        evd["PS2_VeryStrong"],
    ]

    cls = classfyv2(PVS1, PS, PM, PP, BA1, BS, BP, PVS)

    print(f"Evidences: {','.join([i for i in evd if evd[i] == 1])}")
    print(f"ACMG: {cls}")


test()
