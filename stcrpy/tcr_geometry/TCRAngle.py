import os
import math
import numpy as np

from Bio.SVDSuperimposer import SVDSuperimposer


def superimpose(coordsA, coordsB):
    sup = SVDSuperimposer()
    sup.set(coordsA, coordsB)
    sup.run()
    return sup


class TCRAngle:

    def __init__(self, tcr_type="abTCR"):
        """Unified class for calculating TCR orientation angles for abTCRs and gdTCRs."""

        self.dat_path = os.path.join(os.path.dirname(__file__), "reference_data")

        from stcrpy.tcr_processing.TCRParser import TCRParser  # avoids circular import

        self.sparser = TCRParser(QUIET=True)

        self.chain_configs = {
            "abTCR": {
                "chains": ("A", "B"),
                "consensus": ("consensus_A.pdb", "consensus_B.pdb"),
                "pc": ("pcA.txt", "pcB.txt"),
                "coreset": ("Acoreset.txt", "Bcoreset.txt"),
                "angle_labels": ["BA", "BC1", "AC1", "BC2", "AC2", "dc"],
            },
            "gdTCR": {
                "chains": ("G", "D"),
                "consensus": ("consensus_G.pdb", "consensus_D.pdb"),
                "pc": ("pcA.txt", "pcB.txt"),
                "coreset": ("Acoreset.txt", "Bcoreset.txt"),
                "angle_labels": ["DG", "DC1", "GC1", "DC2", "GC2", "dc"],
            },
        }

        self._init_tcr_type_specific_reference_data(tcr_type)

    def _init_tcr_type_specific_reference_data(self, tcr_type):
        if tcr_type not in self.chain_configs:
            raise ValueError("tcr_type must be 'abTCR' or 'gdTCR'")
        self.tcr_type = tcr_type
        self.cfg = self.chain_configs[tcr_type]
        self._read_consensus()
        self._read_pc()
        self._read_coreset()

    def _normalise(self, v):
        a = np.array(v)
        return a / np.linalg.norm(a)

    def _read_consensus(self):
        c1, c2 = self.cfg["chains"]
        f1, f2 = self.cfg["consensus"]
        self.consensus_1 = self.sparser.get_tcr_structure(
            c1, os.path.join(self.dat_path, f1)
        )
        self.consensus_2 = self.sparser.get_tcr_structure(
            c2, os.path.join(self.dat_path, f2)
        )

        self.consensus_1_atoms = sorted(
            list(self.consensus_1.get_atoms()), key=lambda x: x.parent.id[1]
        )
        self.consensus_2_atoms = sorted(
            list(self.consensus_2.get_atoms()), key=lambda x: x.parent.id[1]
        )

    def _read_pc(self):
        f1, f2 = self.cfg["pc"]
        self.pos1 = [
            list(map(float, x.split())) for x in open(os.path.join(self.dat_path, f1))
        ]
        self.pos2 = [
            list(map(float, x.split())) for x in open(os.path.join(self.dat_path, f2))
        ]
        self.c1 = [
            6 * 0.5 * self.pos1[0][i] - 2 * 0.5 * self.pos1[1][i] + self.pos1[2][i]
            for i in range(3)
        ]
        self.c2 = [
            -10 * 0.5 * self.pos2[0][i] + 1 * 0.5 * self.pos2[1][i] + self.pos2[2][i]
            for i in range(3)
        ]
        self.p1a = [self.c1[i] + self.pos1[0][i] for i in range(3)]
        self.p1b = [self.c1[i] + self.pos1[1][i] for i in range(3)]
        self.p2a = [self.c2[i] + self.pos2[0][i] for i in range(3)]
        self.p2b = [self.c2[i] + self.pos2[1][i] for i in range(3)]

    def _read_coreset(self):
        f1, f2 = self.cfg["coreset"]
        self.coreset1 = [
            int(l.strip()[1:]) for l in open(os.path.join(self.dat_path, f1))
        ]
        self.coreset2 = [
            int(l.strip()[1:]) for l in open(os.path.join(self.dat_path, f2))
        ]

    def calculate_angles(self, tcr):
        if self.tcr_type not in str(type(tcr)):
            self.tcr_type = str(type(tcr)).split(".")[-1].split("'>")[0]
            self._init_tcr_type_specific_reference_data(self.tcr_type)

        c1, c2 = self.cfg["chains"]
        tcr_chain_1_coreset_atoms = [
            res["CA"]
            for res in tcr[tcr.get_domain_assignment()["V" + c1]].get_residues()
            if res.id[1] in self.coreset1
        ]
        tcr_chain_2_coreset_atoms = [
            res["CA"]
            for res in tcr[tcr.get_domain_assignment()["V" + c2]].get_residues()
            if res.id[1] in self.coreset2
        ]

        tcr_chain_1_coreset_atoms = sorted(
            tcr_chain_1_coreset_atoms, key=lambda x: x.parent.id[1]
        )
        tcr_chain_2_coreset_atoms = sorted(
            tcr_chain_2_coreset_atoms, key=lambda x: x.parent.id[1]
        )

        rot1, tran1 = superimpose(
            np.asarray([a.coord for a in tcr_chain_1_coreset_atoms]),
            np.asarray([a.coord for a in self.consensus_1_atoms]),
        ).get_rotran()
        rot2, tran2 = superimpose(
            np.asarray([a.coord for a in tcr_chain_2_coreset_atoms]),
            np.asarray([a.coord for a in self.consensus_2_atoms]),
        ).get_rotran()

        points1 = [
            np.dot(self.c1, rot1) + tran1,
            np.dot(self.p1a, rot1) + tran1,
            np.dot(self.p1b, rot1) + tran1,
        ]
        points2 = [
            np.dot(self.c2, rot2) + tran2,
            np.dot(self.p2a, rot2) + tran2,
            np.dot(self.p2b, rot2) + tran2,
        ]

        C = self._normalise([points2[0][i] - points1[0][i] for i in range(3)])
        Cminus = [-x for x in C]

        v1a = self._normalise([points1[1][i] - points1[0][i] for i in range(3)])
        v1b = self._normalise([points1[2][i] - points1[0][i] for i in range(3)])
        v2a = self._normalise([points2[1][i] - points2[0][i] for i in range(3)])
        v2b = self._normalise([points2[2][i] - points2[0][i] for i in range(3)])

        dc = np.linalg.norm(points2[0] - points1[0])

        n_x = np.cross(v1a, C)
        n_y = np.cross(C, n_x)
        tmp1 = self._normalise([0, np.dot(v1a, n_x), np.dot(v1a, n_y)])
        tmp2 = self._normalise([0, np.dot(v2a, n_x), np.dot(v2a, n_y)])

        angle = math.degrees(math.acos(np.dot(tmp1, tmp2)))
        if np.dot(np.cross(tmp1, tmp2), [1, 0, 0]) < 0:
            angle = -angle

        results = [
            angle,
            math.degrees(math.acos(np.dot(v1a, C))),
            math.degrees(math.acos(np.dot(v2a, Cminus))),
            math.degrees(math.acos(np.dot(v1b, C))),
            math.degrees(math.acos(np.dot(v2b, Cminus))),
            dc,
        ]

        return dict(zip(self.cfg["angle_labels"], results))
