import argparse
import pandas as pd
import os
import ast
import sys
import re

from cyvcf2 import VCF
import time
import subprocess

Strength = {
    "Normal": "Normal",
    "VeryStrong": "VeryStrong",
    "Strong": "Strong",
    "Moderate": "Moderate",
    "Supporting": "Supporting",
    "Unmet": "Unmet",
}

CLS = {
    "PAT": "Pathogenic",
    "LP": "Likely pathogenic",
    "VUS": "Uncertain significance",
    "LB": "Likely benign",
    "BEN": "Benign",
}

R_CLS = {
    "Pathogenic": "PAT",
    "Likely pathogenic": "LP",
    "Uncertain significance": "VUS",
    "Likely benign": "LB",
    "Benign": "BEN",
}

# Global variables
clinvar_vcf = None

nas_string = "."

# Global metrics variables
metrics_PAT = {
    "LP": 0,
    "VUS": 0,
    "LB": 0,
    "BEN": 0,
}

metrics_LP = {
    "PAT": 0,
    "VUS": 0,
    "LB": 0,
    "BEN": 0,
}

metrics_VUS = {
    "PAT": 0,
    "LP": 0,
    "LB": 0,
    "BEN": 0,
}

metrics_LB = {
    "PAT": 0,
    "LP": 0,
    "VUS": 0,
    "BEN": 0,
}

metrics_BEN = {
    "PAT": 0,
    "LP": 0,
    "VUS": 0,
    "LB": 0,
}

Special_Evidence_Mapping = {
    "PM3_VeryStrong": {
        "evidence": "PM3",
        "strength": Strength["VeryStrong"],
    },
    "PS2_VeryStrong": {
        "evidence": "PS2",
        "strength": Strength["VeryStrong"],
    },
    "PVS1_Strong": {
        "evidence": "PVS1",
        "strength": Strength["Strong"],
    },
    "PP3_Strong": {
        "evidence": "PP3",
        "strength": Strength["Strong"],
    },
    "PVS1_Moderate": {
        "evidence": "PVS1",
        "strength": Strength["Moderate"],
    },
    "PP3_Moderate": {
        "evidence": "PP3",
        "strength": Strength["Moderate"],
    },
    "PVS1_Supporting": {
        "evidence": "PVS1",
        "strength": Strength["Supporting"],
    },
    "PM2_Supporting": {
        "evidence": "PM2",
        "strength": Strength["Supporting"],
    },
    "BP4_Strong": {
        "evidence": "BP4",
        "strength": Strength["Strong"],
    },
}

# Evidence list
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

Evidence_Types = {
    f"{Strength['VeryStrong']}": ["PVS1", "BA1"],
    f"{Strength['Strong']}": [
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
    f"{Strength['Moderate']}": [
        "PM1",
        "PM2",
        "PM3",
        "PM4",
        "PM5",
        "PM6",
    ],
    f"{Strength['Supporting']}": [
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

Benign_Evidence_Types = {
    f"{Strength['VeryStrong']}": ["BA1"],
    f"{Strength['Strong']}": [
        "BS1",
        "BS2",
        "BS3",
        "BS4",
    ],
    f"{Strength['Supporting']}": [
        "BP1",
        "BP2",
        "BP3",
        "BP4",
        "BP5",
        "BP6",
        "BP7",
    ],
}


def initial_evidences():
    return {
        **{ev: Strength["Unmet"] for ev in Pathogenic_Evidences},
        **{ev: Strength["Unmet"] for ev in Benign_Evidences},
    }


def parse_evidences(evidence_list):
    """
    Convert evidence list like ["PM2", "PVS1_Strong"] to initial_evidences() format.

    Args:
        evidence_list: List of evidence codes, optionally with strength suffixes
                      e.g. ["PM2", "PVS1_Strong", "PP1_Supporting"]

    Returns:
        Dict with evidence codes as keys and strength values
    """
    evidences = initial_evidences()

    for evidence in evidence_list:
        if "_" in evidence:
            # Handle evidence with strength suffix (e.g., "PVS1_Strong")
            code, strength_suffix = evidence.split("_", 1)
            if code in evidences:
                if strength_suffix in Strength.values():
                    evidences[code] = strength_suffix
        else:
            # Handle plain evidence code (e.g., "PM2")
            if evidence in evidences:
                evidences[evidence] = Strength["Normal"]

    return evidences


def evidences_to_str(evidences):
    """
    Convert evidence dictionary to list of evidence codes.
    Reverses the parse_evidences function.

    Args:
        evidences: Dict with evidence codes as keys and strength values
                  e.g. {"PM2": "Normal", "PVS1": "Strong", "PP1": "Unmet"}

    Returns:
        List of evidence codes with strength suffixes for non-default strengths
        e.g. ["PM2", "PVS1_Strong"]
    """
    evidence_list = []

    for evidence, strength in evidences.items():
        if strength == Strength["Unmet"]:
            continue  # Skip unmet evidences

        if strength == Strength["Normal"]:
            # Use plain evidence code for normal strength
            evidence_list.append(evidence)
        else:
            # Use evidence code with strength suffix for other strengths
            evidence_list.append(f"{evidence}_{strength}")

    return ",".join(evidence_list) if len(evidence_list) > 0 else nas_string


def get_intervar_data(intervar):
    """
    Example can be:
    InterVar: Benign PVS1=0 PS=[0, 0, 0, 0, 0] PM=[0, 0, 0, 0, 0, 0, 0] PP=[0, 0, 0, 0, 0, 0] BA1=1 BS=[1, 0, 0, 0, 0] BP=[0, 0, 0, 0, 0, 0, 0, 0]
    Or can be
    InterVar: Benign PVS1=0 PS=[0, 0, 0, 0, 0] PM=[0, 0, 0, 0, 0, 0, 0] PP=[0, 0, 0, 0, 0, 0] BA1=1 BS=[1, 0, 0, 0, 0] BP=[0, 0, 0, 0, 0, 0, 0, 0] PVS=[0, 0]
    Extract original InterVar evidences and return list of active evidence codes
    """
    evidences = initial_evidences()

    if intervar == nas_string:
        return nas_string, evidences

    classification = intervar.split("InterVar: ")[1].split(" PVS1")[0]

    info_list = intervar.split("InterVar: ")[1].split("=")

    PVS1 = int(info_list[1].split(" PS")[0])
    PS = ast.literal_eval(info_list[2].split(" PM")[0])
    PM = ast.literal_eval(info_list[3].split(" PP")[0])
    PP = ast.literal_eval(info_list[4].split(" BA1")[0])
    BA1 = int(info_list[5].split(" BS")[0])
    BS = ast.literal_eval(info_list[6].split(" BP")[0])
    BP = ast.literal_eval(info_list[7].split(" PVS")[0])
    PVS = ast.literal_eval(info_list[8]) if len(info_list) >= 9 else [0, 0]

    # Evidence code mappings
    evidence_mappings = [
        (PVS1 == 1, ["PVS1"]),
        (True, [("PM3_VeryStrong", "PS2_VeryStrong"), PVS]),
        (True, [("PS1", "PS2", "PS3", "PS4", "PVS1_Strong", "PP3_Strong"), PS]),
        (
            True,
            [
                (
                    "PM1",
                    "PM2",
                    "PM3",
                    "PM4",
                    "PM5",
                    "PM6",
                    "PVS1_Moderate",
                    "PP3_Moderate",
                ),
                PM,
            ],
        ),
        (
            True,
            [
                (
                    "PP1",
                    "PP2",
                    "PP3",
                    "PP4",
                    "PP5",
                    "PVS1_Supporting",
                    "PM2_Supporting",
                ),
                PP,
            ],
        ),
        (BA1 == 1, ["BA1"]),
        (True, [("BS1", "BS2", "BS3", "BS4", "BP4_Strong"), BS]),
        (True, [("BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7"), BP]),
    ]

    for condition, mapping in evidence_mappings:
        if not condition:
            continue

        if len(mapping) == 1:  # Single evidence
            if isinstance(mapping[0], list):
                for ev in mapping[0]:
                    evidences[ev] = Strength["Normal"]
            else:
                evidences[mapping[0]] = Strength["Normal"]
        else:  # Code-value pair
            codes, values = mapping
            for i, value in enumerate(values):
                if value == 1 and i < len(codes):
                    evidence_code = codes[i]
                    if evidence_code in Special_Evidence_Mapping:
                        mapped = Special_Evidence_Mapping[evidence_code]
                        evidences[mapped["evidence"]] = mapped["strength"]
                    else:
                        evidences[evidence_code] = Strength["Normal"]

    return classification, evidences


def modify_evidences_based_on_autoPVS1(current_pvs1, autopvs1_strength):
    strength_mapping = {
        Strength["VeryStrong"]: Strength["VeryStrong"],
        Strength["Strong"]: Strength["Strong"],
        Strength["Moderate"]: Strength["Moderate"],
        Strength["Supporting"]: Strength["Supporting"],
        Strength["Unmet"]: Strength["Unmet"],
        nas_string: Strength["Unmet"],
    }

    new_pvs1 = (
        current_pvs1
        if autopvs1_strength not in strength_mapping
        else strength_mapping[autopvs1_strength]
    )
    adjusted = 1 if new_pvs1 != current_pvs1 else 0

    return new_pvs1, adjusted


def check_PM1(chr_pos_ref_alt):
    # Extract chrom, pos, ref, alt from key_data
    key_data = chr_pos_ref_alt.split("_")
    chrom = key_data[0]
    pos = int(key_data[1])

    # Loop through clinvar regions
    start_pos = pos - 25
    end_pos = pos + 25
    region = f"{chrom}:{start_pos}-{end_pos}"

    count = 0

    # Fetch variants within the specified range
    for variant in clinvar_vcf(region):
        try:
            clnsig = variant.INFO["CLNSIG"]
            if "Pathogenic" in clnsig or "Likely_pathogenic" in clnsig:
                count += 1
        except KeyError:
            # No CLNSIG found
            pass
    return Strength["Normal"] if count >= 4 else Strength["Unmet"], count


def check_PP3(REVEL, SpliceAI):
    REVEL = float(REVEL)
    SpliceAI = float(SpliceAI) if SpliceAI != nas_string else -1000
    if (REVEL >= 0.644 and REVEL < 0.773) or (SpliceAI > 0.2):
        return Strength["Supporting"]
    if REVEL >= 0.932:
        return Strength["Strong"]
    if REVEL >= 0.773 and REVEL < 0.932:
        return Strength["Moderate"]
    return Strength["Unmet"]


def check_BP4(REVEL):
    REVEL = float(REVEL)
    if REVEL > 0 and REVEL <= 0.003:
        return Strength["VeryStrong"]
    if REVEL > 0.003 and REVEL <= 0.016:
        return Strength["Strong"]
    if REVEL > 0.016 and REVEL <= 0.290:
        return Strength["Supporting"]
    return Strength["Unmet"]


def get_evidences_by_strength(evidences, strength, type="pathogenic"):
    possible_list = Pathogenic_Evidences if type == "pathogenic" else Benign_Evidences
    return [
        ev
        for ev in evidences.keys()
        if ev in possible_list
        and (
            evidences[ev] == strength
            or (ev in Evidence_Types[strength] and evidences[ev] == Strength["Normal"])
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

    PVS1 = 1 if evidences["PVS1"] == Strength["Normal"] else 0
    BA1 = 1 if evidences["BA1"] == Strength["Normal"] else 0
    BP7 = 1 if evidences["BP7"] == Strength["Normal"] else 0
    BS1 = 1 if evidences["BS1"] == Strength["Normal"] else 0

    PVS = get_evidences_by_strength(evidences, Strength["VeryStrong"])
    PS = get_evidences_by_strength(evidences, Strength["Strong"])
    PM = get_evidences_by_strength(evidences, Strength["Moderate"])
    PP = get_evidences_by_strength(evidences, Strength["Supporting"])

    BS = get_evidences_by_strength(evidences, Strength["Strong"], "benign")
    BP = get_evidences_by_strength(evidences, Strength["Supporting"], "benign")

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
            if evidences["PM2"] == Strength["Supporting"]:
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


def get_priority(cls):
    CLS_Priority = {
        "Pathogenic": 0,
        "Likely pathogenic": 1,
        "Uncertain significance": 2,
        "Likely benign": 3,
        "Benign": 4,
    }
    return 6 if cls == nas_string else CLS_Priority[cls]


def check_truncating_variant(coding_effects, var_locations, cnomen):
    r"""
    Splicing Variant meet one of the following conditions:
        1. Variant with var_locations contain: "splicing"
        2. cnomen with pattern r'^c\.\\d+[\\+\\-]([3-9]|10)[ATCG].*' or pattern r'^c\.\\d+[\\+\\-][12][ATCG].*'
    Variant with coding_effects contain one of the following:
        1. "Frameshift"
        2. "Start loss"
        3. "Start retained"
        4. "Stop gained"
        5. "Stop lost"
        6. "Stop retained"
    """
    # Splicing variant
    if (
        "splicing" in var_locations
        or re.match(r"^c\.\d+[\+\-]([3-9]|10)[ATCG].*", cnomen)
        or re.match(r"^c\.\d+[\+\-][12][ATCG].*", cnomen)
    ):
        return True

    # Coding Effect conditions
    valid_coding_effects = [
        "frameshift",
        "start loss",
        "start retained",
        "stop gained",
        "stop lost",
        "stop retained",
    ]

    # Convert arrays to sets and check for intersection
    if set(valid_coding_effects) & set(coding_effects):
        return True

    return False


def check_missense_variant(coding_effects):
    if "missense" in coding_effects:
        return True
    return False


def check_inframe_indel_variant(coding_effects):
    if "inframe insertion" in coding_effects or "inframe deletion" in coding_effects:
        return True
    return False


def check_synonymous_variant(coding_effects):
    if "synonymous" in coding_effects:
        return True
    return False


def do_restrict_evidences(
    is_truncating,
    is_missense,
    is_inframe_indel,
    is_synonymous,
    oevidences,
):
    evidences = oevidences.copy()
    match_list = []
    """
          Truncating  Missense    Inframe Indel   Synonymous
    * PS1     X           /             X             X
      PS2     /           /             /             /
      PS3     /           /             /             /
      PS4     /           /             /             /
    * PM1     X           /             /             X
      PM2     /           /             /             /
      PM3     /           /             /             /
    * PM4     X           X             /             X
    * PM5     X           /             X             X
      PM6     /           /             /             /
      PP1     /           /             /             /
    * PP2     X           /             X             X
    * PP3     X           /             /             /
      PP4     /           /             /             /
    * PP5     X           X             X             X

      BA1     /           /             /             /
      BS1     /           /             /             /
      BS2     /           /             /             /
      BS3     /           /             /             /
      BS4     /           /             /             /
    * BP1     X           /             X             X
      BP2     /           /             /             /
    * BP3     X           X             /             X
    * BP4     X           /             /             /
      BP5     /           /             /             /
    * BP6     X           X             X             X
    * BP7     X           X             X             /
    """

    # PS1 & PM5 & PP2 & BP1
    if not is_missense and (is_truncating or is_inframe_indel or is_synonymous):
        # PS1
        evidences["PS1"] = Strength["Unmet"]
        # PM5
        evidences["PM5"] = Strength["Unmet"]
        # PP2
        evidences["PP2"] = Strength["Unmet"]
        # BP1
        evidences["BP1"] = Strength["Unmet"]
        match_list.extend(["PS1", "PM5", "PP2", "BP1"])

    # PM1
    if not is_missense and not is_inframe_indel and (is_truncating or is_synonymous):
        evidences["PM1"] = Strength["Unmet"]
        match_list.append("PM1")

    # PM4
    if not is_inframe_indel and (is_truncating or is_missense or is_synonymous):
        # PM4
        evidences["PM4"] = Strength["Unmet"]
        # BP3
        evidences["BP3"] = Strength["Unmet"]
        match_list.extend(["PM4", "BP3"])

    # PM5 same condition with PS1

    # PP2 same condition with PS1

    # PP3
    if is_truncating and not is_missense and not is_inframe_indel and not is_synonymous:
        # PP3
        evidences["PP3"] = Strength["Unmet"]
        # BP4
        evidences["BP4"] = Strength["Unmet"]
        match_list.extend(["PP3", "BP4"])

    # PP5
    if is_truncating and is_missense and is_inframe_indel and is_synonymous:
        # PP5
        evidences["PP5"] = Strength["Unmet"]
        # BP6
        evidences["BP6"] = Strength["Unmet"]
        match_list.extend(["PP5", "BP6"])

    # BP1 same condition with PS1

    # BP3 same condition with PM4

    # BP4 same condition with PP3

    # BP6 same condition with PP5

    # BP7
    if not is_synonymous and (is_truncating or is_missense or is_inframe_indel):
        evidences["BP7"] = Strength["Unmet"]
        match_list.append("BP7")

    return evidences, match_list


def restrict_evidences(variant, oevidences):
    vcoding_effect = variant["coding_effect"]
    vvar_location = variant["var_location"]
    cnomen = variant["c_nomen"]

    coding_effects = [e.lower() for e in vcoding_effect.split(", ")]
    var_locations = [e.lower() for e in vvar_location.split(", ")]

    # Check Truncating variant
    is_truncating = check_truncating_variant(coding_effects, var_locations, cnomen)

    # Check Missense variant
    is_missense = check_missense_variant(coding_effects)

    # Check Inframe Indel variant
    is_inframe_indel = check_inframe_indel_variant(coding_effects)

    # Check Synonymous variant
    is_synonymous = check_synonymous_variant(coding_effects)

    # Restrict evidences based on variant type
    evidences, match_list = do_restrict_evidences(
        is_truncating, is_missense, is_inframe_indel, is_synonymous, oevidences
    )

    return evidences, match_list


def format_list(list):
    return ", ".join([str(e) for e in list])


def get_InterVar_str(cls, PVS1, PS, PM, PP, BA1, BS, BP, PVS):
    return f"InterVar: {cls} PVS1={PVS1} PS=[{format_list(PS)}] PM=[{format_list(PM)}] PP=[{format_list(PP)}] BA1={BA1} BS=[{format_list(BS)}] BP=[{format_list(BP)}] PVS=[{format_list(PVS)}]"


def modify_denovo_intervar_info(row):
    """
    Add PS2 evidence to every denovo variants
    """

    # Evidences list:
    PVS1_standalone = "PVS1"
    PVS_list = ["PM3_VeryStrong", "PS2_VeryStrong"]
    PS_list = ["PS1", "PS2", "PS3", "PS4", "PVS1_Strong", "PP3_Strong"]
    PM_list = [
        "PM1",
        "PM2",
        "PM3",
        "PM4",
        "PM5",
        "PM6",
        "PVS1_Moderate",
        "PP3_Moderate",
    ]
    PP_list = ["PP1", "PP2", "PP3", "PP4", "PP5", "PVS1_Supporting", "PM2_Supporting"]
    BA1_standalone = "BA1"
    BS_list = ["BS1", "BS2", "BS3", "BS4", "BP4_Strong"]
    BP_list = ["BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7"]

    InterVar = row["InterVar"]

    InterVar_str_modified = nas_string

    ACMG_evidences = nas_string
    ACMG_cls = nas_string
    ACMG_priority = nas_string

    if "ACMG_evidences" in row:
        evidences_strength = parse_evidences(
            (
                row["ACMG_evidences"].split(",")
                if row["ACMG_evidences"] != nas_string
                else []
            )
        )
        evidences_strength["PS2"] = Strength["Normal"]
        ACMG_cls = classify(evidences=evidences_strength)
        ACMG_priority = get_priority(ACMG_cls)
        ACMG_evidences = evidences_to_str(evidences_strength)
    else:

        PVS1 = nas_string
        PVS = nas_string
        PS = nas_string
        PM = nas_string
        PP = nas_string
        BA1 = nas_string
        BS = nas_string
        BP = nas_string

        if InterVar != nas_string:
            InterVar_info_list = InterVar.split("InterVar: ")[1].split("=")

            PVS1 = int(InterVar_info_list[1].split(" PS")[0])
            PS = ast.literal_eval(InterVar_info_list[2].split(" PM")[0])
            PM = ast.literal_eval(InterVar_info_list[3].split(" PP")[0])
            PP = ast.literal_eval(InterVar_info_list[4].split(" BA1")[0])
            BA1 = int(InterVar_info_list[5].split(" BS")[0])
            BS = ast.literal_eval(InterVar_info_list[6].split(" BP")[0])
            BP = ast.literal_eval(InterVar_info_list[7].split(" PVS")[0])
            if len(InterVar_info_list) >= 9:
                PVS = ast.literal_eval(InterVar_info_list[8])
            else:
                PVS = [0 for _ in PVS_list]

            # Add PS2 evidence
            PS[1] = 1

            ACMG_cls = classify(PVS1, PS, PM, PP, BA1, BS, BP, PVS)
            InterVar_str_modified = get_InterVar_str(
                ACMG_cls, PVS1, PS, PM, PP, BA1, BS, BP, PVS
            )
            ACMG_priority = get_priority(ACMG_cls)

        else:
            PVS1 = 0
            PS = [0 for _ in PS_list]
            PM = [0 for _ in PM_list]
            PP = [0 for _ in PP_list]
            BA1 = 0
            BS = [0 for _ in BS_list]
            BP = [0 for _ in BP_list]
            PVS = [0 for _ in PVS_list]

            # Add PS2 evidence
            PS[1] = 1

            ACMG_cls = classify(PVS1, PS, PM, PP, BA1, BS, BP, PVS)
            InterVar_str_modified = get_InterVar_str(
                ACMG_cls, PVS1, PS, PM, PP, BA1, BS, BP, PVS
            )
            ACMG_priority = get_priority(ACMG_cls)

    result = {
        "key": row["chr_pos_ref_alt_gene"],
        "InterVar": InterVar_str_modified,
        "ACNG_evidences": ACMG_evidences,
        "ACMG_classification": ACMG_cls,
        "ACMG_priority": ACMG_priority,
    }

    return result


def modify_intervar_info(row):
    InterVar = row["InterVar"]
    strength = row["AutoPVS1_strength"]
    REVEL = row["REVEL_score"]
    SpliceAI = row["SpliceAI_max_score"]
    key = row["chr_pos_ref_alt"]

    AutoPVS1_adjusted = 0
    PP3_modified = 0
    BP4_modified = 0

    InterVar_cls, InterVar_evidences = get_intervar_data(InterVar)
    InterVar_priority = get_priority(InterVar_cls)

    ACMG_modified_evidences = InterVar_evidences.copy()

    ##### Start Modify intervar evidences #####
    # 1. PVS1 modified by AutoPVS1 strength
    ACMG_modified_evidences["PVS1"], AutoPVS1_adjusted = (
        modify_evidences_based_on_autoPVS1(InterVar_evidences["PVS1"], strength)
    )

    # 2. Change all PM2 to PM2_supporting
    if InterVar_evidences["PM2"] != Strength["Unmet"]:
        ACMG_modified_evidences["PM2"] = Strength["Supporting"]

    # 3. PM1 Eg. For A mutation, atleast 4 missense variants that have P or LP or P/LP on either side of A (including A itself)
    ACMG_modified_evidences["PM1"], _ = check_PM1(key)

    # 4. PP3, BP4
    # Check PP3
    PP3_strength = check_PP3(REVEL, SpliceAI)
    PP3_modified = 1 if PP3_strength != Strength["Unmet"] else 0
    ACMG_modified_evidences["PP3"] = PP3_strength

    # Check BP4
    BP4_strength = check_BP4(REVEL)
    BP4_modified = 1 if BP4_strength != Strength["Unmet"] else 0
    # Trigger BA1 if BP4 is VeryStrong
    if BP4_strength == Strength["VeryStrong"]:
        ACMG_modified_evidences["BA1"] = Strength["Normal"]
        ACMG_modified_evidences["BP4"] = Strength["Unmet"]
    else:
        ACMG_modified_evidences["BP4"] = BP4_strength

    ACMG_modified_cls = classify(ACMG_modified_evidences)
    ACMG_modified_priority = get_priority(ACMG_modified_cls)

    ##### End Modify intervar evidences #####

    ##### Start restrict evidences #####
    ACMG_final_evidences = ACMG_modified_evidences.copy()

    ACMG_final_evidences, restrict_list = restrict_evidences(
        row, ACMG_modified_evidences
    )
    ACMG_final_cls = classify(ACMG_final_evidences)
    ACMG_final_priority = get_priority(ACMG_final_cls)
    ##### End restrict evidences #####

    result = {
        "key": row["chr_pos_ref_alt"] + "_" + row["gene"],
        "AutoPVS1_adjusted": AutoPVS1_adjusted,
        # Original InterVar
        "InterVar_cls": InterVar_cls,
        "InterVar_evidences": evidences_to_str(InterVar_evidences),
        "InterVar_priority": InterVar_priority,
        # Modified InterVar
        "ACMG_modified_cls": ACMG_modified_cls,
        "ACMG_modified_evidences": evidences_to_str(ACMG_modified_evidences),
        "ACMG_modified_priority": ACMG_modified_priority,
        "PP3_modified": PP3_modified,
        "BP4_modified": BP4_modified,
        # Modified & Restricted InterVar
        "ACMG_final_cls": ACMG_final_cls,
        "ACMG_final_evidences": evidences_to_str(ACMG_final_evidences),
        "ACMG_final_priority": ACMG_final_priority,
        "restrict_list": ",".join(restrict_list) if len(restrict_list) else nas_string,
    }

    return result


def write_logs(results):
    """
    Write comprehensive logs for classification analysis.

    Args:
        results: List of result dictionaries from modify_intervar_info
    """
    from collections import defaultdict, Counter

    # Initialize counters and matrices
    intervar_to_modified_matrix = defaultdict(lambda: defaultdict(int))
    intervar_to_final_matrix = defaultdict(lambda: defaultdict(int))

    pp3_modified_count = 0
    bp4_modified_count = 0
    pp3_records = []
    bp4_records = []
    restricted_records = []

    cls_counts = {
        "InterVar_cls": Counter(),
        "ACMG_modified_cls": Counter(),
        "ACMG_final_cls": Counter(),
    }

    # Process each result
    for result in results:
        intervar_cls = result["InterVar_cls"]
        modified_cls = result["ACMG_modified_cls"]
        final_cls = result["ACMG_final_cls"]

        # Update classification matrices
        intervar_to_modified_matrix[intervar_cls][modified_cls] += 1
        intervar_to_final_matrix[intervar_cls][final_cls] += 1

        # Count PP3 and BP4 modifications
        if result["PP3_modified"] == 1:
            pp3_modified_count += 1
            pp3_records.append(result)

        if result["BP4_modified"] == 1:
            bp4_modified_count += 1
            bp4_records.append(result)

        # Process restricted evidences
        final_evidences = (
            result["ACMG_final_evidences"].split(",")
            if result["ACMG_final_evidences"] != nas_string
            else []
        )
        restrict_list = (
            result["restrict_list"].split(",")
            if result["restrict_list"] != nas_string
            else []
        )

        # Find restricted evidences that exist in final evidences
        restricted_evidences = [ev for ev in final_evidences if ev in restrict_list]

        restricted_records.append(
            {
                "key": result["key"],
                "ACMG_final_evidences": result["ACMG_final_evidences"],
                "restrict_list": result["restrict_list"],
                "restricted_list": (
                    ",".join(restricted_evidences)
                    if restricted_evidences
                    else nas_string
                ),
            }
        )

        # Count classifications
        cls_counts["InterVar_cls"][intervar_cls] += 1
        cls_counts["ACMG_modified_cls"][modified_cls] += 1
        cls_counts["ACMG_final_cls"][final_cls] += 1

    # Write classification matrices in readable format
    with open("cls_matrix.log", "w") as f:
        f.write("Classification Change Matrix: InterVar_cls -> ACMG_modified_cls\n")
        f.write("=" * 80 + "\n")

        # Get all classifications - ensure all 5 classifications are included
        all_cls = [
            "Pathogenic",
            "Likely pathogenic",
            "Uncertain significance",
            "Likely benign",
            "Benign",
        ]

        # Write header with abbreviated column names
        col_abbrev = {
            "Pathogenic": "PAT",
            "Likely pathogenic": "LP",
            "Uncertain significance": "VUS",
            "Likely benign": "LB",
            "Benign": "BEN",
        }
        f.write(f"{'From/To':<25}")
        for cls in all_cls:
            f.write(f"{col_abbrev[cls]:<8}")
        f.write("Total\n")
        f.write("-" * 80 + "\n")

        # Write matrix with row totals (only show differences)
        for from_cls in all_cls:
            row_total = sum(
                intervar_to_modified_matrix[from_cls][to_cls]
                for to_cls in all_cls
                if from_cls != to_cls
            )
            f.write(f"{from_cls:<25}")
            for to_cls in all_cls:
                count = intervar_to_modified_matrix[from_cls][to_cls]
                # Only show count if it's a change (from_cls != to_cls) or if no change but count > 0
                if from_cls != to_cls:
                    f.write(f"{count:<8}")
                else:
                    f.write(f"{'.':<8}")  # Use dot for diagonal (no change)
            f.write(f"{row_total}\n")

        # Write column totals (only show changes)
        f.write("-" * 80 + "\n")
        f.write(f"{'Total':<25}")
        changes_total = 0
        for to_cls in all_cls:
            col_total = sum(
                intervar_to_modified_matrix[from_cls][to_cls]
                for from_cls in all_cls
                if from_cls != to_cls
            )
            f.write(f"{col_total:<8}")
            changes_total += col_total
        f.write(f"{changes_total}\n")

        f.write("\n\n")
        f.write("Classification Change Matrix: InterVar_cls -> ACMG_final_cls\n")
        f.write("=" * 80 + "\n")

        # Write header for second matrix
        f.write(f"{'From/To':<25}")
        for cls in all_cls:
            f.write(f"{col_abbrev[cls]:<8}")
        f.write("Total\n")
        f.write("-" * 80 + "\n")

        # Write second matrix with row totals (only show differences)
        for from_cls in all_cls:
            row_total = sum(
                intervar_to_final_matrix[from_cls][to_cls]
                for to_cls in all_cls
                if from_cls != to_cls
            )
            f.write(f"{from_cls:<25}")
            for to_cls in all_cls:
                count = intervar_to_final_matrix[from_cls][to_cls]
                # Only show count if it's a change (from_cls != to_cls) or if no change but count > 0
                if from_cls != to_cls:
                    f.write(f"{count:<8}")
                else:
                    f.write(f"{'.':<8}")  # Use dot for diagonal (no change)
            f.write(f"{row_total}\n")

        # Write column totals for second matrix (only show changes)
        f.write("-" * 80 + "\n")
        f.write(f"{'Total':<25}")
        changes_total = 0
        for to_cls in all_cls:
            col_total = sum(
                intervar_to_final_matrix[from_cls][to_cls]
                for from_cls in all_cls
                if from_cls != to_cls
            )
            f.write(f"{col_total:<8}")
            changes_total += col_total
        f.write(f"{changes_total}\n")

    # Write PP3 log
    with open("PP3.log", "w") as f:
        f.write(f"PP3 Modified Records: {pp3_modified_count} total\n")
        f.write("=" * 50 + "\n")
        f.write("key\tInterVar_cls\tACMG_modified_cls\tACMG_final_cls\n")
        for record in pp3_records:
            f.write(
                f"{record['key']}\t{record['InterVar_cls']}\t{record['ACMG_modified_cls']}\t{record['ACMG_final_cls']}\n"
            )

    # Write BP4 log
    with open("BP4.log", "w") as f:
        f.write(f"BP4 Modified Records: {bp4_modified_count} total\n")
        f.write("=" * 50 + "\n")
        f.write("key\tInterVar_cls\tACMG_modified_cls\tACMG_final_cls\n")
        for record in bp4_records:
            f.write(
                f"{record['key']}\t{record['InterVar_cls']}\t{record['ACMG_modified_cls']}\t{record['ACMG_final_cls']}\n"
            )

    # Write restricted log
    with open("restricted.log", "w") as f:
        # Calculate total records with restrictions applied
        restricted_count = sum(
            1
            for record in restricted_records
            if record["restricted_list"] != nas_string
        )

        f.write("Restricted Evidences Analysis\n")
        f.write("=" * 50 + "\n")
        f.write(
            f"Total records with restricted evidences applied: {restricted_count}\n"
        )
        f.write("=" * 50 + "\n")
        f.write("key\tACMG_final_evidences\trestrict_list\trestricted_list\n")
        for record in restricted_records:
            f.write(
                f"{record['key']}\t{record['ACMG_final_evidences']}\t{record['restrict_list']}\t{record['restricted_list']}\n"
            )

    # Write total classification counts
    with open("cls_total.log", "w") as f:
        f.write("Total Classification Counts\n")
        f.write("=" * 30 + "\n")

        for cls_type, counts in cls_counts.items():
            f.write(f"\n{cls_type}:\n")
            f.write("-" * 20 + "\n")
            for classification, count in sorted(counts.items()):
                f.write(f"{classification}: {count}\n")

    print(f"Logs written:")
    print(f"- Classification matrices: cls_matrix.log")
    print(f"- PP3 modifications: PP3.log ({pp3_modified_count} records)")
    print(f"- BP4 modifications: BP4.log ({bp4_modified_count} records)")
    print(f"- Restricted evidences: restricted.log")
    print(f"- Total classification counts: cls_total.log")


def run_denovo(input, output):
    start_time = time.time()

    # Load input file
    df = pd.read_csv(input, sep="\t")

    # Apply the function to the DataFrame
    parsed_data = df.apply(modify_denovo_intervar_info, axis=1)

    # Convert the list of dictionaries to a DataFrame
    parsed_df = pd.DataFrame(parsed_data.tolist())

    # Export result
    parsed_df.to_csv(output, sep="\t", index=False)

    print(f"Whole process took: {time.time() - start_time} secs")


def run(input, output, clinvar):
    global clinvar_vcf

    start_time = time.time()

    # calculate_metrics(output=output)
    # sys.exit()

    # Load Clinvar VCF
    clinvar_vcf = VCF(clinvar)

    # Load input file
    df = pd.read_csv(input, sep="\t")

    # Apply the function to the DataFrame
    parsed_data = df.apply(modify_intervar_info, axis=1)

    # Write comprehensive logs
    parsed_data_list = parsed_data.tolist()
    write_logs(parsed_data_list)

    # Convert the list of dictionaries to a DataFrame
    parsed_df = pd.DataFrame(parsed_data_list)

    # Export result
    parsed_df.to_csv(output, sep="\t", index=False)

    print("")
    print(f"Whole process took: {time.time() - start_time} secs")


if __name__ == "__main__":
    # python3 modify_InterVar.py -i example.tsv -o out.tsv -c clinvar_missense.vcf.gz

    parser = argparse.ArgumentParser(
        description="""Modify InterVar evidences & classification result V2"""
    )
    parser.add_argument("-i", "--input", type=str, help="Input tsv file", required=True)
    parser.add_argument(
        "-c",
        "--clinvar",
        dest="clinvar",
        type=str,
        help="Clinvar VCF file",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        type=str,
        help="Output file with modified result",
        required=True,
    )
    parser.add_argument(
        "-d", "--denovo", dest="denovo", type=bool, help="Denovo analysis"
    )
    args = parser.parse_args()

    if args.denovo:
        print("Running denovo modification")
        run_denovo(input=args.input, output=args.output)
    else:
        print("Running modification")
        run(input=args.input, output=args.output, clinvar=args.clinvar)
