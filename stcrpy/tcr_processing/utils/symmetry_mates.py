import warnings
import tempfile
import os

from Bio.PDB import PDBParser, PDBIO
from ..TCRParser import TCRParser


def get_symmetry_mates(filename):
    try:
        from pymol import cmd
    except ModuleNotFoundError:
        warnings.warn(
            "Pymol not installed - please install pymol to enabe symmetry mate TCR parsing. "
        )
        return []

    cmd.load(filename)
    obj_name = cmd.get_object_list()[0]
    cmd.symexp(obj_name, obj_name, obj_name, cutoff=20.0, quiet=0)
    if (
        len(cmd.get_object_list()) == 1
    ):  # No symmetry mates found, likely becuase file did not contain symmetry information
        cmd.delete(obj_name)
        cmd.fetch(
            filename.split("/")[-1].split(".")[0]
        )  # this will try to retrieve the file from the pdb directly, will not work if the filename is not the pdb code
    if len(cmd.get_object_list()) == 0:
        warnings.warn(f"No symmetry mates found for {filename}.")
        return
    obj_name = cmd.get_object_list()[0]
    cmd.symexp(obj_name, obj_name, obj_name, cutoff=10.0, quiet=0)
    tcr_symmetry_mates = []

    tcp = TCRParser()
    pdp = PDBParser(QUIET=True)
    pdbio = PDBIO()
    with tempfile.TemporaryDirectory() as tmpdir:
        for i, obj in enumerate(cmd.get_object_list()):
            fn = os.path.join(tmpdir, f"{obj_name}_symmetry_mate_{i}.pdb")
            cmd.save(fn, obj)
            symmetry_mate = pdp.get_structure("tmp", fn)
            if i == 0:
                chain_ids = generate_chain_id_list(
                    len(list(symmetry_mate.get_chains()) * len(cmd.get_object_list()))
                )
                for c in symmetry_mate.get_chains():
                    chain_ids.remove(
                        c.id
                    )  # remove all chain ids of the original structure
            if i > 0:  # Skip the original structure
                # rename chain ids, this cannot be done directly to TCR structure without breaking the TCR and MHC chain assignments.
                for chain in reversed(list(symmetry_mate.get_chains())):
                    symmetry_mate[0].detach_child(chain.id)
                    new_id = chain_ids.pop(0)
                    chain.id = new_id
                    symmetry_mate[0].add(chain)
                pdbio.set_structure(symmetry_mate)
                pdbio.save(fn)

                symmetry_mate = tcp.get_tcr_structure(
                    f"{obj_name}_symmetry_{i}", fn, include_symmetry_mates=False
                )
                tcr_symmetry_mates.append(symmetry_mate)

    # clean up the pymol cmd space
    for obj in cmd.get_object_list():
        cmd.delete(obj)
    del cmd
    return tcr_symmetry_mates


def generate_chain_id_list(N):
    """Generates a set of chain ids starting from A, B, C, ..., Z, AA, AB, ..., AZ, BA

    Args:
        N (int): The number of chain IDs to generate.

    Returns:
        set: A set of generated chain IDs.
    """
    chain_ids = []
    for i in range(N):
        ch1 = chr(65 + (i // 26) - 1) if i >= 26 else ""
        ch2 = chr(65 + (i % 26))
        chain_ids.append(ch1 + ch2)
    import string

    chain_ids = string.ascii_uppercase + string.ascii_lowercase + string.digits
    return list(chain_ids)[:N]
