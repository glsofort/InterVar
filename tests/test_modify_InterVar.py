import unittest
import modify_InterVar as m


class TestingFunc:
    # def test_check_format_1(self, cases, func):
    #     for case in cases:
    #         result = func(
    #             coding_effects=case["coding_effects"],
    #             var_locations=case["var_locations"],
    #             cnomen=case["cnomen"],
    #         )
    #         self.assertEqual(result, case["result"], f"Test ID: {case['id']}")
    def restrict_evidences(str_list, evd_list, evds):
        tmp_evd_list = evd_list[:]
        for evd in evds:
            if evd not in str_list:
                continue
            index = str_list.index(evd)
            tmp_evd_list[index] = 0

        return tmp_evd_list


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
                "PS": [1, 1, 0, 0, 0, 0],
                "PM": [0, 1, 0, 0, 0, 0, 0, 0],
                "PP": [0, 0, 1, 0, 0, 0, 0],
                "BP": [0, 0, 1, 0, 0, 0, 0],
                "result_PS": [],
                "result_PM": [],
                "result_PP": [],
                "result_BP": [],
                "result_evidences": ["PS1", "PM5", "PP2", "BP1"],
            },
            {
                "id": 2,
                "is_truncating": False,
                "is_missense": False,
                "is_inframe_indel": True,
                "is_synonymous": True,
                "PS": [1, 1, 0, 0, 0, 0],
                "PM": [0, 1, 0, 0, 0, 0, 0, 0],
                "PP": [0, 0, 1, 0, 0, 0, 0],
                "BP": [0, 0, 1, 0, 0, 0, 0],
                "result_PS": [],
                "result_PM": [],
                "result_PP": [],
                "result_BP": [],
                "result_evidences": ["PS1", "PM5", "PP2", "BP1"],
            },
            {
                "id": 3,
                "is_truncating": False,
                "is_missense": False,
                "is_inframe_indel": True,
                "is_synonymous": False,
                "PS": [1, 1, 0, 0, 0, 0],
                "PM": [0, 1, 0, 0, 0, 0, 0, 0],
                "PP": [0, 0, 1, 0, 0, 0, 0],
                "BP": [0, 0, 1, 0, 0, 0, 0],
                "result_PS": [],
                "result_PM": [],
                "result_PP": [],
                "result_BP": [],
                "result_evidences": ["PS1", "PM5", "PP2", "BP1", "BP7"],
            },
            {
                "id": 4,
                "is_truncating": True,
                "is_missense": False,
                "is_inframe_indel": False,
                "is_synonymous": False,
                "PS": [1, 1, 0, 0, 0, 0],
                "PM": [0, 1, 0, 0, 0, 0, 0, 0],
                "PP": [0, 0, 1, 0, 0, 0, 0],
                "BP": [0, 0, 1, 0, 0, 0, 0],
                "result_PS": [],
                "result_PM": [],
                "result_PP": [],
                "result_BP": [],
                "result_evidences": [
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
            },
            {
                "id": 5,
                "is_truncating": True,
                "is_missense": True,
                "is_inframe_indel": False,
                "is_synonymous": False,
                "PS": [1, 1, 0, 0, 0, 0],
                "PM": [0, 1, 0, 0, 0, 0, 0, 0],
                "PP": [0, 0, 1, 0, 0, 0, 0],
                "BP": [0, 0, 1, 0, 0, 0, 0],
                "result_PS": [],
                "result_PM": [],
                "result_PP": [],
                "result_BP": [],
                "result_evidences": ["PM4", "BP3", "BP7"],
            },
        ]

        for case in cases:
            # Update test result list
            case["result_PS"] = TestingFunc.restrict_evidences(
                m.PS_list, case["PS"], case["result_evidences"]
            )
            case["result_PM"] = TestingFunc.restrict_evidences(
                m.PM_list, case["PM"], case["result_evidences"]
            )
            case["result_PP"] = TestingFunc.restrict_evidences(
                m.PP_list, case["PP"], case["result_evidences"]
            )
            case["result_BP"] = TestingFunc.restrict_evidences(
                m.BP_list, case["BP"], case["result_evidences"]
            )

            # Test
            result_PS, result_PM, result_PP, result_BP, result_evidences = (
                m.do_restrict_evidences(
                    is_truncating=case["is_truncating"],
                    is_missense=case["is_missense"],
                    is_inframe_indel=case["is_inframe_indel"],
                    is_synonymous=case["is_synonymous"],
                    PS=case["PS"],
                    PM=case["PM"],
                    PP=case["PP"],
                    BP=case["BP"],
                )
            )

            # Assertions
            self.assertCountEqual(
                result_evidences,
                case["result_evidences"],
                f"Test ID: {case['id']}. Case: All.",
            )
            self.assertCountEqual(
                result_PS, case["result_PS"], f"Test ID: {case['id']}. Case: PS."
            )
            self.assertCountEqual(
                result_PM, case["result_PM"], f"Test ID: {case['id']}. Case: PM."
            )
            self.assertCountEqual(
                result_PP, case["result_PP"], f"Test ID: {case['id']}. Case: PP."
            )
            self.assertCountEqual(
                result_BP, case["result_BP"], f"Test ID: {case['id']}. Case: BP."
            )

    def test_restrict_evidences(self):
        pass


if __name__ == "__main__":
    unittest.main()
