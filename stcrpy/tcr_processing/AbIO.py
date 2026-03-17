from Bio import PDB
from Bio.PDB.PDBIO import PDBIO

from .Ab import Ab


class AbIO(PDBIO):
    def __init__(self):
        super().__init__()
        self.io = PDBIO()

    def save(
        self,
        ab: Ab,
        save_as: str = None,
        ab_only: bool = False,
        format: str = "pdb",
    ):
        assert isinstance(ab, Ab), f"{ab} must be type Ab"

        structure_to_save = PDB.Model.Model(0)

        for chain in ab.get_chains():
            chain.serial_num = 0
            structure_to_save.add(chain)

        if not ab_only:
            for chain in ab.get_antigen():
                chain.serial_num = 0
                structure_to_save.add(chain)

        self.io.set_structure(structure_to_save)

        if not save_as:
            parent_id = ab.parent.parent.id if ab.parent and ab.parent.parent else ab.id
            if not ab_only:
                save_as = f"{parent_id}_{ab.id}.{format}"
            else:
                save_as = f"{parent_id}_{ab.id}_Ab_only.{format}"

        self.io.save(save_as)
