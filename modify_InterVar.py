import argparse
import pandas as pd
import os
import ast
import sys
from cyvcf2 import VCF
import time
import subprocess

# Constant
evidence_keys = ["PVS1", "PS", "PM", "PP", "BA1", "BS", "BP"]
nas_string = "."
Strength = {
    "VeryStrong": "VeryStrong",
    "Strong": "Strong",
    "Moderate": "Moderate",
    "Supporting": "Supporting",
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


def sum_of_list(list):
    sum = 0
    for i in list:
        sum = sum + i
    return sum


def classify(PVS1, PS, PM, PP, BA1, BS, BP):
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

    if (PVS1 == 1 or PS_sum > 0 or PM_sum > 0 or PP_sum > 0) and (
        BA1 == 1 or BS_sum > 0 or BP_sum > 0
    ):
        return BPS[4]  # Uncertain significance x

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
        return BPS[0]  # Pathogenic x

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


def modify_evidences_based_on_autoPVS1(PVS1, PS, PM, PP, AutoPVS1_strength):
    adjusted = 0
    # PVS1
    if AutoPVS1_strength == Strength["VeryStrong"]:
        if PVS1 != 1:
            adjusted = 1
        PVS1 = 1

    # Add PVS1_strong
    if AutoPVS1_strength == Strength["Strong"]:
        PVS1 = 0
        PS[4] = 1
        adjusted = 1

    # Add PVS1_moderate
    if AutoPVS1_strength == Strength["Moderate"]:
        PVS1 = 0
        PM[6] = 1
        adjusted = 1

    # Add PVS1_supporting
    if AutoPVS1_strength == Strength["Supporting"]:
        PVS1 = 0
        PP[5] = 1
        adjusted = 1

    return PVS1, PS, PM, PP, adjusted


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
    return 1 if count >= 4 else 0, count


def check_PP3(REVEL, SpliceAI):
    REVEL = float(REVEL)
    SpliceAI = float(SpliceAI) if SpliceAI != nas_string else -1000
    if REVEL >= 0.644 and REVEL < 0.773 and SpliceAI > 0.2:
        return Strength["Supporting"]
    if REVEL >= 0.932:
        return Strength["Strong"]
    if REVEL >= 0.773 and REVEL < 0.932:
        return Strength["Moderate"]
    return nas_string


def check_BP4(REVEL):
    REVEL = float(REVEL)
    if REVEL > 0 and REVEL <= 0.003:
        return Strength["VeryStrong"]
    if REVEL > 0.003 and REVEL <= 0.016:
        return Strength["Strong"]
    if REVEL > 0.016 and REVEL <= 0.290:
        return Strength["Supporting"]
    return nas_string


def format_list(list):
    return ", ".join([str(e) for e in list])


def get_evidences(intervar):
    if intervar == nas_string:
        return nas_string

    info_list = intervar.split("InterVar: ")[1].split("=")

    PVS1 = int(info_list[1].split(" PS")[0])
    PS = ast.literal_eval(info_list[2].split(" PM")[0])
    PM = ast.literal_eval(info_list[3].split(" PP")[0])
    PP = ast.literal_eval(info_list[4].split(" BA1")[0])
    BA1 = int(info_list[5].split(" BS")[0])
    BS = ast.literal_eval(info_list[6].split(" BP")[0])
    BP = ast.literal_eval(info_list[7])

    evidences = []

    if PVS1 == 1:
        evidences.append("PVS1")

    for i in range(len(PS)):
        if PS[i] == 1:
            if i == 4:
                evidences.append("PVS1_strong")
            elif i == 5:
                evidences.append("PP3_strong")
            else:
                evidences.append(f"PS{i+1}")

    for i in range(len(PM)):
        if PM[i] == 1:
            if i == 6:
                evidences.append("PVS1_moderate")
            elif i == 7:
                evidences.append("PP3_moderate")
            else:
                evidences.append(f"PM{i+1}")

    for i in range(len(PP)):
        if PP[i] == 1:
            if i == 5:
                evidences.append("PVS1_supporting")
            if i == 6:
                evidences.append("PM2_supporting")
            else:
                evidences.append(f"PP{i+1}")

    if BA1 == 1:
        evidences.append("BA1")

    for i in range(len(BS)):
        if BS[i] == 1:
            if i == 4:
                evidences.append(f"BP4_strong")
            else:
                evidences.append(f"BS{i+1}")

    for i in range(len(BP)):
        if BP[i] == 1:
            evidences.append(f"BP{i+1}")

    return ",".join(evidences)


def get_priority(cls):
    CLS_Priority = {
        "Pathogenic": 0,
        "Likely pathogenic": 1,
        "Uncertain significance": 2,
        "Likely benign": 3,
        "Benign": 4,
    }
    return CLS_Priority[cls]


def modify_intervar_info(row):
    InterVar = row["InterVar"]
    strength = row["AutoPVS1_strength"]
    REVEL = row["REVEL_score"]
    SpliceAI = row["SpliceAI_max_score"]
    key = row["chr_pos_ref_alt"]

    AutoPVS1_adjusted = 0
    PP3_modified = 0
    BP4_modified = 0

    InterVar_cls = nas_string
    InterVar_cls_modified = nas_string

    InterVar_modified = nas_string

    InterVar_priority = 6
    InterVar_priority_modified = 6

    InterVar_evidences = get_evidences(InterVar)
    InterVar_evidences_modified = nas_string

    PVS1 = nas_string
    PS = nas_string
    PM = nas_string
    PP = nas_string
    BA1 = nas_string
    BS = nas_string
    BP = nas_string

    # Extract InterVar evidences
    # InterVar: Benign PVS1=0 PS=[0, 0, 0, 0, 0] PM=[0, 0, 0, 0, 0, 0, 0] PP=[0, 0, 0, 0, 0, 0] BA1=1 BS=[1, 0, 0, 0, 0] BP=[0, 0, 0, 0, 0, 0, 0, 0]
    if InterVar != nas_string:
        InterVar_info_list = InterVar.split("InterVar: ")[1].split("=")
        InterVar_cls = InterVar.split("InterVar: ")[1].split(" PVS1")[0]
        InterVar_priority = get_priority(InterVar_cls)

        PVS1 = int(InterVar_info_list[1].split(" PS")[0])
        PS = ast.literal_eval(InterVar_info_list[2].split(" PM")[0])
        PM = ast.literal_eval(InterVar_info_list[3].split(" PP")[0])
        PP = ast.literal_eval(InterVar_info_list[4].split(" BA1")[0])
        BA1 = int(InterVar_info_list[5].split(" BS")[0])
        BS = ast.literal_eval(InterVar_info_list[6].split(" BP")[0])
        BP = ast.literal_eval(InterVar_info_list[7])

        # Append additional evidences in each group
        # Pathogenic
        # PVS1_strong
        PS[4] = 0  # index: 4
        # PP3_strong
        PS.append(0)  # index: 5
        # PVS1_moderate
        PM[6] = 0  # index: 6
        # PP3_moderate
        PM.append(0)  # index: 7
        # PVS1_supporting
        PP[5] = 0  # index: 5
        # PM2_supporting
        PP.append(0)  # index: 6

        # Benign
        # BP4_strong
        BS[4] = 0  # index: 4

        # Modify intervar evidences
        # 1. PVS1 modified by AutoPVS1 strength
        PVS1, PS, PM, PP, AutoPVS1_adjusted = modify_evidences_based_on_autoPVS1(
            PVS1, PS, PM, PP, strength
        )

        # 2. Change all PM2 to PM2_supporting
        PM2 = PM[1]
        if PM2 == 1:
            # Remove PM2
            PM[1] = 0
            # Add PM2_supporting
            PP[6] = 1

        # 3. PM1 Eg. For A mutation, atleast 4 missense variants that have P or LP or P/LP on either side of A (including A itself)
        PM1, _ = check_PM1(key)
        PM[0] = PM1

        # 4. PP3, BP4
        # Check PP3
        PP3_strength = check_PP3(REVEL, SpliceAI)
        if PP3_strength != nas_string:
            PP3_modified = 1
            PP[2] = 0
        if PP3_strength == Strength["Strong"]:
            PS[5] = 1
        elif PP3_strength == Strength["Moderate"]:
            PM[7] = 1
        elif PP3_strength == Strength["Supporting"]:
            PP[2] = 1
        # Check BP4
        BP4_strength = check_BP4(REVEL)
        if BP4_strength != nas_string:
            BP4_modified = 1
            BP[3] = 0
        # Trigger BA1
        if BP4_strength == Strength["VeryStrong"]:
            BA1 = 1
        elif BP4_strength == Strength["Strong"]:
            BS[4] = 1
        elif BP4_strength == Strength["Supporting"]:
            BP[3] = 1

        # Classification
        InterVar_cls_modified = classify(PVS1, PS, PM, PP, BA1, BS, BP)

        # InterVar modification
        # InterVar: Benign PVS1=0 PS=[0, 0, 0, 0, 0, 0] PM=[0, 0, 0, 0, 0, 0, 0, 0] PP=[0, 0, 0, 0, 0, 0, 0, 0] BA1=1 BS=[1, 0, 0, 0, 0] BP=[0, 0, 0, 0, 0, 0, 0, 0]
        InterVar_modified = f"InterVar: {InterVar_cls_modified} PVS1={PVS1} PS=[{format_list(PS)}] PM=[{format_list(PM)}] PP=[{format_list(PP)}] BA1={BA1} BS=[{format_list(BS)}] BP=[{format_list(BP)}]"

        InterVar_evidences_modified = get_evidences(InterVar_modified)

        InterVar_priority_modified = get_priority(InterVar_cls_modified)

    result = {
        "key": row["chr_pos_ref_alt"] + "_" + row["gene"],
        "AutoPVS1_adjusted": AutoPVS1_adjusted,
        "InterVar_cls": InterVar_cls,
        "InterVar_cls_modified": InterVar_cls_modified,
        "InterVar_evidences": InterVar_evidences,
        "InterVar_evidences_modified": InterVar_evidences_modified,
        "PP3_modified": PP3_modified,
        "BP4_modified": BP4_modified,
        "InterVar_priority": InterVar_priority,
        "InterVar_priority_modified": InterVar_priority_modified,
        "InterVar": InterVar,
        "InterVar_modified": InterVar_modified,
    }
    return result


def run_cmd(cmd):
    # Run a simple command, e.g., 'ls' or 'dir' depending on your OS
    result = subprocess.run([cmd], shell=True, capture_output=True, text=True)
    # Print the output of the command
    return result.stdout.strip()


def calculate_metrics(output):
    current_dir = os.getcwd()
    result = os.path.join(current_dir, output)

    # Classifications
    for k, v in R_CLS.items():
        print(f"{v} to others:")
        oth = {}
        for k2, v2 in R_CLS.items():
            if k2 != k:
                cmd = (
                    'awk -F"\t" \'{if($3 != $4 && $3 == "'
                    + k
                    + '" && $4 == "'
                    + k2
                    + '") { print $1"\t"$3"\t"$4 }}\' '
                    + result
                    + " | wc -l"
                )
                res = run_cmd(cmd)
                oth[v2] = int(res)
        print(oth)

    # PP3 & BP4
    cmd = 'awk -F"\t" \'{if($7 == 1) { print $5"\t"$6 }}\' ' + result + " | wc -l"
    res = run_cmd(cmd)
    print(f"Total modified PP3: {res}")

    cmd = 'awk -F"\t" \'{if($8 == 1) { print $5"\t"$6 }}\' ' + result + " | wc -l"
    res = run_cmd(cmd)
    print(f"Total modified BP4: {res}")

    # Add more log file:
    cmd = "awk -F\"\t\" '{if($7 == 1) { print $0 }}' " + result + " > PP3.log"
    run_cmd(cmd)
    cmd = "awk -F\"\t\" '{if($8 == 1) { print $0 }}' " + result + " > BP4.log"
    run_cmd(cmd)


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

    # Convert the list of dictionaries to a DataFrame
    parsed_df = pd.DataFrame(parsed_data.tolist())

    # Export result
    parsed_df.to_csv(output, sep="\t", index=False)

    # Metrics
    print("Metrics:")
    calculate_metrics(output=output)

    print("")
    print(f"Whole process took: {time.time() - start_time} secs")


def test():
    # awk -F"\t" '{if($3 != $4) { print $1"\t"$3"\t"$4 }}' out.tsv | less

    # less out.tsv | awk -F"\t" '{if(index($6, "PM1") != 0){print}}'  | wc -l
    # less out.tsv | awk -F"\t" '{if(index($5, "PM1") != 0){print}}'  | wc -l

    # Number of VUS replaced by others
    # awk -F"\t" '{if($3 != $4 && $3 == "Uncertain significance") { print $1"\t"$3"\t"$4 }}' out.tsv
    # awk -F"\t" '{if($3 != $4 && $3 == "Pathogenic") { print $1"\t"$3"\t"$4 }}' out.tsv

    # Number of others replaced by VUS

    pass


if __name__ == "__main__":
    # python3 modify_InterVar.py -i example.tsv -o out.tsv -c clinvar_missense.vcf.gz

    parser = argparse.ArgumentParser(
        description="""Modify InterVar evidences & classification result"""
    )
    parser.add_argument("-i", "--input", type=str, help="Input tsv file", required=True)
    parser.add_argument(
        "-c",
        "--clinvar",
        dest="clinvar",
        type=str,
        help="Clinvar VCF file",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        type=str,
        help="Output file with modified result",
        required=True,
    )
    args = parser.parse_args()

    run(input=args.input, output=args.output, clinvar=args.clinvar)
