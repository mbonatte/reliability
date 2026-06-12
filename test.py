import copy
import unittest
import warnings
from math import pi

import numpy as np
import secan as sa

from reliability import Reliability


def merge_dict(variables, new_variable):
    variables["names"].append(new_variable["names"])
    variables["dists"].append(new_variable["dists"])
    variables["bounds"].append(new_variable["bounds"])
    return variables


class TestSections(unittest.TestCase):
    def setUp(self):
        variables = self._build_variables()
        beam = self._build_beam()

        self.beam_reliability = Reliability(copy.deepcopy(beam))
        self.beam_reliability.variables = copy.deepcopy(variables)
        self.beam_reliability.inverted = False

    def _build_variables(self):
        variables = {"names": [], "dists": [], "bounds": []}

        fck = 18e6
        cov = 0.2
        bias_factor = 1.30
        rusch = 0.85
        fcm = bias_factor * fck
        std = cov * fcm
        variables = merge_dict(
            variables,
            {"names": "fc", "dists": "norm", "bounds": [rusch * fcm, rusch * std]},
        )

        fyk = 500e6
        cov = 0.04
        bias_factor = 1.22
        fym = bias_factor * fyk
        std = cov * fym
        variables = merge_dict(
            variables,
            {"names": "fy", "dists": "norm", "bounds": [fym, std]},
        )

        variables = merge_dict(
            variables,
            {"names": "As", "dists": "norm", "bounds": [1, 0.02]},
        )
        variables = merge_dict(
            variables,
            {"names": "cover_bottom", "dists": "norm", "bounds": [0, 5e-3]},
        )
        variables = merge_dict(
            variables,
            {"names": "cover_top", "dists": "norm", "bounds": [10e-3, 10e-3]},
        )

        b_web = 0.4
        h_web = 1.6
        h_flange = 0.3

        variables = merge_dict(
            variables,
            {
                "names": "b_web",
                "dists": "norm",
                "bounds": [b_web, (4 + 6 * b_web) * 1e-3],
            },
        )
        variables = merge_dict(
            variables,
            {
                "names": "h_web",
                "dists": "norm",
                "bounds": [h_web, (4 + 6 * h_web) * 1e-3],
            },
        )
        variables = merge_dict(
            variables,
            {
                "names": "h_flange",
                "dists": "norm",
                "bounds": [h_flange, (4 + 6 * h_flange) * 1e-3],
            },
        )

        variables["num_vars"] = len(variables["names"])
        return variables

    def _build_beam(self):
        fck = 18e6
        fyk = 500e6
        gamma_c = 1.4
        gamma_s = 1.15
        rusch = 0.85

        b_web = 0.4
        h_web = 1.6
        b_flange = 3.5
        h_flange = 0.3

        concrete_d = sa.material.Concrete(fc=rusch * fck / gamma_c)
        steel_d = sa.material.SteelIdeal(
            young=210e9,
            fy=fyk / gamma_s,
            ultimate_strain=10e-3,
        )

        web = sa.geometry.RectSection(
            width=b_web,
            height=h_web,
            material=concrete_d,
            center=(0, h_web / 2),
            n_discret=1,
        )
        flange = sa.geometry.RectSection(
            width=b_flange,
            height=h_flange,
            material=concrete_d,
            center=(0, h_web + h_flange / 2),
            n_discret=30,
        )

        beam = sa.Section(section=[web, flange])

        stirrup_diameter = 9.53e-3
        concrete_cover = 3e-2
        spacing_between_layers = 2e-2
        rebar_diameter = 19.05e-3

        h = concrete_cover + stirrup_diameter + rebar_diameter / 2
        for _ in range(4):
            beam.addLineRebar(
                rebar_diameter,
                steel_d,
                spacing=0.03,
                position=[[-0.17, h], [-0.08, h]],
            )
            beam.addLineRebar(
                rebar_diameter,
                steel_d,
                spacing=0.03,
                position=[[0.17, h], [0.08, h]],
            )
            h += rebar_diameter + spacing_between_layers

        beam.addSingleRebar(rebar_diameter, steel_d, position=[-0.17, h])
        beam.addSingleRebar(rebar_diameter, steel_d, position=[0.17, h])

        height = beam.get_section_height()
        h = height - (concrete_cover + stirrup_diameter + rebar_diameter / 2)
        beam.addLineRebar(
            rebar_diameter,
            steel_d,
            spacing=0.03,
            position=[[-0.17, h], [-0.08, h]],
        )
        beam.addLineRebar(
            rebar_diameter,
            steel_d,
            spacing=0.03,
            position=[[0.17, h], [0.08, h]],
        )

        return beam

    def test_variables_match_reference_case(self):
        expected_variables = {
            "names": [
                "fc",
                "fy",
                "As",
                "cover_bottom",
                "cover_top",
                "b_web",
                "h_web",
                "h_flange",
            ],
            "dists": ["norm", "norm", "norm", "norm", "norm", "norm", "norm", "norm"],
            "bounds": [
                [19890000.0, 3978000.0],
                [610000000.0, 24400000.0],
                [1, 0.02],
                [0, 0.005],
                [0.01, 0.01],
                [0.4, 0.0064],
                [1.6, 0.013600000000000001],
                [0.3, 0.0058],
            ],
            "num_vars": 8,
        }

        self.assertEqual(self.beam_reliability.variables, expected_variables)

    def test_reliability_stores_reference_rebar_properties(self):
        expected_area = pi * 19.05e-3**2 / 4

        self.assertEqual(len(self.beam_reliability.rebar_area), 42)
        self.assertEqual(len(self.beam_reliability.rebar_pos), 42)
        for area in self.beam_reliability.rebar_area:
            self.assertAlmostEqual(area, expected_area)

    def test_update_rebar_area_can_be_reset_to_original(self):
        original_areas = list(self.beam_reliability.rebar_area)

        self.beam_reliability.update_beam("As", 1.10)
        updated_areas = [
            section.area
            for section in self.beam_reliability.beam.section
            if isinstance(section, sa.geometry.Rebar)
        ]

        for updated_area, original_area in zip(updated_areas, original_areas):
            self.assertAlmostEqual(updated_area, original_area * 1.10)

        self.beam_reliability.update_beam("As", 1.10)
        updated_again_areas = [
            section.area
            for section in self.beam_reliability.beam.section
            if isinstance(section, sa.geometry.Rebar)
        ]

        for updated_area, original_area in zip(updated_again_areas, original_areas):
            self.assertAlmostEqual(updated_area, original_area * 1.10)

        self.beam_reliability.set_original()
        reset_areas = [
            section.area
            for section in self.beam_reliability.beam.section
            if isinstance(section, sa.geometry.Rebar)
        ]

        for reset_area, original_area in zip(reset_areas, original_areas):
            self.assertAlmostEqual(reset_area, original_area)

    def test_get_resistance_moment_returns_reference_value(self):
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="invalid value encountered in scalar divide",
                category=RuntimeWarning,
            )
            moment = self.beam_reliability.get_resistance_moment(n_points=10)

        self.assertAlmostEqual(moment, 7235374.981867483)

    def test_sensitivity_bk_returns_reference_values(self):
        self.beam_reliability.n_points = 5

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning)
            sensitivity = self.beam_reliability.sensitivity_bk()

        expected = np.array(
            [
                8.57140908,
                100.0,
                50.54845617,
                0.0,
                2.02204229,
                4.60588899e-10,
                17.14459431,
                7.30622907,
            ]
        )

        self.assertEqual(sensitivity.shape, (self.beam_reliability.variables["num_vars"],))
        np.testing.assert_allclose(sensitivity, expected, rtol=1e-6)

    def test_sensitivity_sobol_returns_seeded_reference_values(self):
        self.beam_reliability.n_points = 5

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning)
            s1, st = self.beam_reliability.sensitivity_Sobol(
                n=2,
                seed=123,
                calc_second_order=False,
            )

        expected_s1 = np.array(
            [
                1.39747542e-01,
                2.27096560e00,
                2.84898360e-01,
                -9.98496354e-03,
                6.40418949e-03,
                -3.28776844e-09,
                -1.45288235e-01,
                -2.20906514e-01,
            ]
        )
        expected_st = np.array(
            [
                6.69100168e-03,
                1.34797002e00,
                2.23200742e-02,
                8.67411032e-03,
                1.03420339e-04,
                4.62271224e-18,
                5.84128175e-03,
                1.38404588e-02,
            ]
        )

        self.assertEqual(s1.shape, (self.beam_reliability.variables["num_vars"],))
        self.assertEqual(st.shape, (self.beam_reliability.variables["num_vars"],))
        np.testing.assert_allclose(s1, expected_s1, rtol=1e-6, atol=1e-9)
        np.testing.assert_allclose(st, expected_st, rtol=1e-6, atol=1e-12)

    def test_monte_carlo_simplified_returns_seeded_reference_values(self):
        self.beam_reliability.n_points = 5

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning)
            moments = self.beam_reliability.monte_carlo_simplified(
                n_LHS=3,
                seed=123,
            )

        expected = np.array(
            [
                10161790.98841443,
                9333964.45228934,
                10733708.82886019,
            ]
        )

        self.assertEqual(moments.shape, (3,))
        np.testing.assert_allclose(moments, expected, rtol=1e-9)


if __name__ == "__main__":
    unittest.main()
