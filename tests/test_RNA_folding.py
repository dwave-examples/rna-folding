# Copyright 2021 D-Wave Systems
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os
import subprocess
import sys
import unittest

import RNA_folding
from dwave.system import LeapHybridCQMSampler

project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class TestSmoke(unittest.TestCase):
    @unittest.skipIf(os.getenv('SKIP_INT_TESTS'), "Skipping integration test.")
    def test_smoke(self):
        """Run RNA_folding.py and check that nothing crashes"""

        demo_file = os.path.join(project_dir, 'RNA_folding.py')
        subprocess.check_output([sys.executable, demo_file])


class TestRNA_folding(unittest.TestCase):
    def test_read_file_to_stem_dict(self):
        """Test ability to read file and create appropriate stem dict"""

        file = os.path.join(project_dir, 'RNA_text_files/simple.txt')

        bond_matrix = RNA_folding.text_to_matrix(file, 2)
        stem_dict = RNA_folding.make_stem_dict(bond_matrix, 3, 2)

        self.assertEqual(stem_dict, {(2, 5, 12, 15): [(2, 4, 13, 15), (2, 5, 12, 15), (3, 5, 12, 14)]})

    def test_build_cqm(self):
        """Test build_CQM creating correct constraints and variables"""

        stem_dict = {
            (1, 3, 13, 15): [(1, 3, 13, 15)],
            (6, 10, 20, 24): [(6, 8, 22, 24), (6, 9, 21, 24), (6, 10, 20, 24), (7, 9, 21, 23), (7, 10, 20, 23),
                              (8, 10, 20, 22)],
            (7, 9, 14, 16): [(7, 9, 14, 16)],
            (13, 15, 23, 25): [(13, 15, 23, 25)]
        }

        cqm = RNA_folding.build_cqm(stem_dict, 3, 0.3)

        self.assertEqual(cqm.variables,
                         [(1, 3, 13, 15), (6, 8, 22, 24), (6, 9, 21, 24), (6, 10, 20, 24), (7, 9, 21, 23),
                          (7, 10, 20, 23),
                          (8, 10, 20, 22), (7, 9, 14, 16), (13, 15, 23, 25), 'Null:(6, 10, 20, 24)']
                         )
        target_linear = {(1, 3, 13, 15): -9.0, (6, 8, 22, 24): -9.0, (6, 9, 21, 24): -16.0, (6, 10, 20, 24): -25.0,
                         (7, 9, 21, 23): -9.0, (7, 10, 20, 23): -16.0, (8, 10, 20, 22): -9.0, (7, 9, 14, 16): -9.0,
                         (13, 15, 23, 25): -9.0, 'Null:(6, 10, 20, 24)': 0.0}
        for v, bias in target_linear.items():
            self.assertEqual(cqm.objective.get_linear(v), bias)

        self.assertEqual(cqm.objective.quadratic,
                         {((6, 8, 22, 24), (1, 3, 13, 15)): 2.6999999999999997,
                          ((6, 9, 21, 24), (1, 3, 13, 15)): 3.5999999999999996, ((6, 10, 20, 24), (1, 3, 13, 15)): 4.5,
                          ((7, 9, 21, 23), (1, 3, 13, 15)): 2.6999999999999997,
                          ((7, 10, 20, 23), (1, 3, 13, 15)): 3.5999999999999996,
                          ((8, 10, 20, 22), (1, 3, 13, 15)): 2.6999999999999997}
                         )

        self.assertEqual(len(cqm.constraints), 15)

    def test_small_case(self):
        """Test solution quality of small case."""

        file = os.path.join(project_dir, 'RNA_text_files/simple_pseudo.txt')

        bond_matrix = RNA_folding.text_to_matrix(file, 2)
        stem_dict = RNA_folding.make_stem_dict(bond_matrix, 3, 2)
        cqm = RNA_folding.build_cqm(stem_dict, 3, 0.3)

        sampler = LeapHybridCQMSampler()
        sample_set = sampler.sample_cqm(cqm)
        stems = RNA_folding.process_cqm_solution(sample_set)

        self.assertEqual(set(stems), {(1, 3, 13, 15), (6, 10, 20, 24)})
