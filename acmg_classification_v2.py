# Constants
nas_string = "."

# Strength
Strength = {
    "Normal": "Normal",
    "VeryStrong": "VeryStrong",
    "Strong": "Strong",
    "Moderate": "Moderate",
    "Supporting": "Supporting",
    "Unmet": "Unmet",
}

# Evidence types
Evidence_Types = {
    "VeryStrong": ["PVS1", "BA1"],
    "Strong": [
        # Pathogenic
        "PS1",
        "PS2",
        "PS3",
        "PS4",
        # Benign
        "BS1",
        "BS2",
        "BS3",
        "BS4",
    ],
    "Moderate": [
        "PM1",
        "PM2",
        "PM3",
        "PM4",
        "PM5",
        "PM6",
    ],
    "Supporting": [
        # Pathogenic
        "PP1",
        "PP2",
        "PP3",
        "PP4",
        "PP5",
        # Benign
        "BP1",
        "BP2",
        "BP3",
        "BP4",
        "BP5",
        "BP6",
        "BP7",
    ],
}

# Pathogenic evidence list
Pathogenic_Evidences = [
    # PVS
    "PVS1",
    # PS
    "PS1",
    "PS2",
    "PS3",
    "PS4",
    # PM
    "PM1",
    "PM2",
    "PM3",
    "PM4",
    "PM5",
    "PM6",
    # PP
    "PP1",
    "PP2",
    "PP3",
    "PP4",
    "PP5",
]

# Benign evidence list
Benign_Evidences = [
    # BA
    "BA1",
    # BS
    "BS1",
    "BS2",
    "BS3",
    "BS4",
    # BP
    "BP1",
    "BP2",
    "BP3",
    "BP4",
    "BP5",
    "BP6",
    "BP7",
]


def get_evidences_list_by_strength(evidences, strength, type="pathogenic"):
    possible_list = Pathogenic_Evidences if type == "pathogenic" else Benign_Evidences
    return [
        ev
        for ev, ev_strength in evidences
        if ev in possible_list
        and (
            ev_strength == strength
            or (ev in Evidence_Types[strength] and ev_strength == nas_string)
        )
    ]


def classify(evidences):
    BPS = [
        "Pathogenic",
        "Likely pathogenic",
        "Benign",
        "Likely benign",
        "Uncertain significance",
    ]

    PVS1 = 1 if "PVS1" in evidences else 0
    BA1 = 1 if "BA1" in evidences else 0
    BP7 = 1 if "BP7" in evidences else 0
    BS1 = 1 if "BS1" in evidences else 0
    PM2_Supporting = 1 if "PM2_Supporting" in evidences else 0

    evidences = [
        (ev.split("_")[0], ev.split("_")[1]) if "_" in ev else (ev, nas_string)
        for ev in evidences
    ]

    PVS = get_evidences_list_by_strength(evidences, Strength["VeryStrong"])
    PS = get_evidences_list_by_strength(evidences, Strength["Strong"])
    PM = get_evidences_list_by_strength(evidences, Strength["Moderate"])
    PP = get_evidences_list_by_strength(evidences, Strength["Supporting"])

    BS = get_evidences_list_by_strength(evidences, Strength["Strong"], "benign")
    BP = get_evidences_list_by_strength(evidences, Strength["Supporting"], "benign")

    PVS_sum = len(PVS)
    PS_sum = len(PS)
    PM_sum = len(PM)
    PP_sum = len(PP)

    BS_sum = len(BS)
    BP_sum = len(BP)

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
            and PVS_sum == 1
            and PS_sum == 0
            and PM_sum == 0
            and PP_sum == 1  # 1 because PM2_Supporting is exist
            and BA1 == 0
            and BS_sum == 0
            and BP_sum == 0
        ):
            if PM2_Supporting == 1:
                return BPS[1]
    except KeyError:
        pass

    if (PVS_sum > 0 or PS_sum > 0 or PM_sum > 0 or PP_sum > 0) and (
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
