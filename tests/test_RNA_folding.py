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
    def test_num_vars(self):
        """Test CQM characteristics for demo instance"""

        file = os.path.join(project_dir, 'RNA_text_files/TMGMV_UPD-PK1.txt')

        bond_matrix = RNA_folding.text_to_matrix(file, 2)
        stem_dict = RNA_folding.make_stem_dict(bond_matrix, 3, 2)
        cqm = RNA_folding.build_cqm(stem_dict, 3, 0.3)

        self.assertEqual(len(cqm.variables), 17)
        self.assertEqual(len(cqm.constraints), 36)

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
