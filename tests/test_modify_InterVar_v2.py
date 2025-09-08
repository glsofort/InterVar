import unittest
import modify_InterVar_v2 as m


class TestingFunc:
    pass


class TestRestrictionFunc(unittest.TestCase):
    def test_check_truncating_variant(self):
        # c.132+2G>A
        # c.1322G>A
        # c.132+10G>A
        cases = [
            {
                "id": 0,
                "coding_effects": ["frameshift"],
                "var_locations": ["exonic"],
                "cnomen": "c.1322G>A",
                "result": True,
            },
            {
                "id": 1,
                "coding_effects": ["incomplete terminal codon"],
                "var_locations": ["exonic"],
                "cnomen": "c.1322G>A",
                "result": False,
            },
            {
                "id": 2,
                "coding_effects": ["incomplete terminal codon", "splice donor"],
                "var_locations": ["splicing", "exonic"],
                "cnomen": "c.1322G>A",
                "result": True,
            },
            {
                "id": 3,
                "coding_effects": ["incomplete terminal codon"],
                "var_locations": ["exonic"],
                "cnomen": "c.132+2G>A",
                "result": True,
            },
            {
                "id": 4,
                "coding_effects": ["incomplete terminal codon"],
                "var_locations": ["exonic"],
                "cnomen": "c.132+10G>A",
                "result": True,
            },
            {
                "id": 5,
                "coding_effects": ["frameshift", "incomplete terminal codon"],
                "var_locations": ["exonic"],
                "cnomen": "c.1322G>A",
                "result": True,
            },
        ]

        for case in cases:
            result = m.check_truncating_variant(
                coding_effects=case["coding_effects"],
                var_locations=case["var_locations"],
                cnomen=case["cnomen"],
            )
            self.assertEqual(result, case["result"], f"Test ID: {case['id']}")

    def test_check_missense_variant(self):
        pass

    def test_check_inframe_indel_variant(self):
        pass

    def test_check_synonymous_variant(self):
        pass

    def test_do_restrict_evidences(self):
        cases = [
            {
                "id": 1,
                "is_truncating": True,
                "is_missense": False,
                "is_inframe_indel": True,
                "is_synonymous": True,
                "test_evidences": ["PS1", "PS2", "PM2", "PP3", "BP3"],
                "result_match_list": ["PS1", "PM5", "PP2", "BP1"],
                "result_evidences": ["PS2", "PM2", "PP3", "BP3"],
            },
            {
                "id": 2,
                "is_truncating": False,
                "is_missense": False,
                "is_inframe_indel": True,
                "is_synonymous": True,
                "test_evidences": ["PS1", "PS2", "PM2", "PP3", "BP3"],
                "result_match_list": ["PS1", "PM5", "PP2", "BP1"],
                "result_evidences": ["PS2", "PM2", "PP3", "BP3"],
            },
            {
                "id": 3,
                "is_truncating": False,
                "is_missense": False,
                "is_inframe_indel": True,
                "is_synonymous": False,
                "test_evidences": ["PS1", "PS2", "PM2", "PP3", "BP3"],
                "result_match_list": ["PS1", "PM5", "PP2", "BP1", "BP7"],
                "result_evidences": ["PS2", "PM2", "PP3", "BP3"],
            },
            {
                "id": 4,
                "is_truncating": True,
                "is_missense": False,
                "is_inframe_indel": False,
                "is_synonymous": False,
                "test_evidences": ["PS1", "PS2", "PM2", "PP3", "BP3"],
                "result_match_list": [
                    "PS1",
                    "PM5",
                    "PP2",
                    "BP1",
                    "BP7",
                    "PP3",
                    "BP4",
                    "PM4",
                    "BP3",
                    "PM1",
                ],
                "result_evidences": ["PS2", "PM2"],
            },
            {
                "id": 5,
                "is_truncating": True,
                "is_missense": True,
                "is_inframe_indel": False,
                "is_synonymous": False,
                "test_evidences": ["PS1", "PS2", "PM2", "PP3", "BP3"],
                "result_match_list": ["PM4", "BP3", "BP7"],
                "result_evidences": [
                    "PS1",
                    "PS2",
                    "PM2",
                    "PP3",
                ],
            },
        ]

        for case in cases:
            # Update test result list
            case["test_evidences"] = m.parse_evidences(case["test_evidences"])
            case["result_evidences"] = m.parse_evidences(case["result_evidences"])

            # Test
            result_evidences, match_list = m.do_restrict_evidences(
                is_truncating=case["is_truncating"],
                is_missense=case["is_missense"],
                is_inframe_indel=case["is_inframe_indel"],
                is_synonymous=case["is_synonymous"],
                oevidences=case["test_evidences"],
            )

            # Assertions
            self.assertCountEqual(
                match_list,
                case["result_match_list"],
                f"Test ID: {case['id']}. Case: Lists should contain same elements.",
            )

            # print("----------------------------------------")
            # print(result_evidences)
            # print(case["result_evidences"])
            # print("----------------------------------------")

            self.assertEqual(
                result_evidences,
                case["result_evidences"],
                f"Test ID: {case['id']}. Evidence dictionaries must match exactly.",
            )

    def test_restrict_evidences(self):
        pass

    def test_classification(self):
        cases = [
            {
                "id": 1,
                "evidences": ["PS1", "PS2", "PM2", "PP3"],
                "result_cls": m.CLS["PAT"],
            },
            {
                "id": 2,
                "evidences": ["PVS1", "PM2_Supporting"],
                "result_cls": m.CLS["LP"],
            },
            {
                "id": 3,
                "evidences": ["PVS1", "PM2_Supporting", "PS1"],
                "result_cls": m.CLS["PAT"],
            },
            {
                "id": 4,
                "evidences": ["PVS1", "PM1_Supporting"],
                "result_cls": m.CLS["VUS"],
            },
            {
                "id": 5,
                "evidences": ["PVS1", "PS1_Supporting"],
                "result_cls": m.CLS["VUS"],
            },
            {
                "id": 6,
                "evidences": ["PVS1", "PS1_Moderate"],
                "result_cls": m.CLS["LP"],
            },
            {
                "id": 7,
                "evidences": ["PVS1_Strong", "PS1_Moderate"],
                "result_cls": m.CLS["LP"],
            },
            {
                "id": 8,
                "evidences": ["PVS1_Strong", "PM2_Supporting"],
                "result_cls": m.CLS["VUS"],
            },
            {
                "id": 9,
                "evidences": ["PVS1_Strong", "PM1_Strong"],
                "result_cls": m.CLS["PAT"],
            },
            {
                "id": 10,
                "evidences": ["PVS1_Strong", "PM2_VeryStrong"],
                "result_cls": m.CLS["PAT"],
            },
            {
                "id": 10,
                "evidences": ["PVS1_Supporting", "PM2_VeryStrong"],
                "result_cls": m.CLS["VUS"],
            },
            {
                "id": 11,
                "evidences": ["PVS1_Moderate", "PM2_VeryStrong"],
                "result_cls": m.CLS["LP"],
            },
            {
                "id": 12,
                "evidences": ["PVS1_Moderate", "PM2_VeryStrong", "BA1"],
                "result_cls": m.CLS["VUS"],
            },
            {
                "id": 13,
                "evidences": ["BA1"],
                "result_cls": m.CLS["BEN"],
            },
            {
                "id": 14,
                "evidences": ["BA1", "BS1"],
                "result_cls": m.CLS["BEN"],
            },
            {
                "id": 15,
                "evidences": ["BA1", "BS1_Supporting"],
                "result_cls": m.CLS["BEN"],
            },
            {
                "id": 16,
                "evidences": ["BS1_Supporting"],
                "result_cls": m.CLS["VUS"],
            },
            {
                "id": 17,
                "evidences": ["BS1"],
                "result_cls": m.CLS["LB"],
            },
            {
                "id": 18,
                "evidences": ["BS1", "BP1", "BP2"],
                "result_cls": m.CLS["LB"],
            },
            {
                "id": 19,
                "evidences": ["BS1", "BP1", "BP2", "BP3"],
                "result_cls": m.CLS["BEN"],
            },
            {
                "id": 20,
                "evidences": ["BS1", "BP7"],
                "result_cls": m.CLS["LB"],
            },
            {
                "id": 21,
                "evidences": ["BS2", "BP1", "BP2"],
                "result_cls": m.CLS["LB"],
            },
            {
                "id": 22,
                "evidences": ["BS2", "BP1", "BP7"],
                "result_cls": m.CLS["BEN"],
            },
        ]

        for case in cases:
            # Update test result list
            case["test_evidences"] = m.parse_evidences(case["evidences"])

            # Test
            result_cls = m.classify(case["test_evidences"])

            if result_cls != case["result_cls"]:
                print(case["test_evidences"])

            # Assertions
            self.assertEqual(
                result_cls,
                case["result_cls"],
                f"Test ID: {case['id']}. Evidence dictionaries must match exactly.",
            )


class TestModifyInterVar(unittest.TestCase):
    def test(self):
        pass


if __name__ == "__main__":
    unittest.main()
