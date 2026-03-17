"""
Created on 3rd April 2024
Adapted from TCR class to Antibody class
"""

import warnings

from Bio import BiopythonWarning

from .Entity import Entity
from .TCRchain import TCRchain  # consider renaming to ABchain if you add one
from .utils.region_definitions import IMGT_VARIABLE_DOMAIN


class ABError(Exception):
    """Error raised when there is an issue with the antibody."""


class Ab(Entity):
    """
    Antibody base class. Inherits from Entity.
    This is a base class for antibody structures, enabling antigen association.
    hlAb is the instantiated subclass of this class.
    """

    def _add_antigen(self, antigen=None):
        """
        Append associated antigen to antibody antigen field.

        Args:
            antigen: Antigen to associate with antibody.
        """
        if antigen not in self.antigen:
            self.antigen.append(antigen)

    def copy(self, copy_siblings=True):
        """
        Return a shallow copy of the AB object.

        If copy_siblings is True, associated antigen objects are also copied.
        If False, antigen references are shared with the original.

        Args:
            copy_siblings: Whether to copy associated antigen objects.

        Returns:
            AB: copied antibody object
        """
        shallow = super().copy()
        if copy_siblings:
            shallow.antigen = [a.copy() for a in self.get_antigen()]
        return shallow

    def get_antigen(self):
        """
        Return a list of antibody-associated antigens.
        """
        return self.antigen

    def is_bound(self):
        """
        True if the antibody is associated with an antigen.

        Returns:
            bool
        """
        return bool(self.get_antigen())

    def get_chains(self):
        """
        Returns generator of antibody chains.
        """
        for c in self:
            yield c

    def get_residues(self):
        """
        Returns generator of antibody residues.
        """
        for c in self.get_chains():
            for r in c:
                yield r

    def get_atoms(self):
        """
        Returns generator of antibody atoms.
        """
        for r in self.get_residues():
            for a in r:
                yield a

    def get_frameworks(self):
        """
        Obtain framework regions from an antibody structure object as generator.
        """
        for f in self.get_fragments():
            if "fw" in f.id:
                yield f

    def get_CDRs(self):
        """
        Obtain complementarity determining regions (CDRs) from an antibody structure object as generator.
        """
        for f in self.get_fragments():
            if "cdr" in f.id:
                yield f

    def get_AB_type(self):
        """
        Get antibody type according to variable region assignments.

        Returns:
            str: antibody type
        """
        if hasattr(self, "ab_type"):
            return self.ab_type
        elif hasattr(self, "VH") and hasattr(self, "VL"):
            self.ab_type = "hlAb"
            return self.ab_type
        return None

    def get_germline_assignments(self):
        """
        Retrieve germline assignments for all antibody chains.

        Returns:
            dict: dict with chain ID as key and germline assignments as value
        """
        return {c.id: c.get_germline_assignments() for c in self.get_chains()}

    def get_germlines_and_antigen(self):
        """
        Get germline assignments for antibody chains together with chain sequences.
        Antigen is included as sequence only and is not assumed to be processed via ANARCI.

        Returns:
            dict: Dictionary of antibody germlines and associated sequences.
        """
        from ..tcr_formats.tcr_formats import get_sequences

        germlines_and_antigen = {}

        try:
            germlines = self.get_germline_assignments()
            for ab_domain, c in self.get_domain_assignment().items():
                germlines_and_antigen[ab_domain] = (
                    germlines[c]["v_gene"][0][1],
                    germlines[c]["j_gene"][0][1],
                )
                germlines_and_antigen[f"{ab_domain}_species"] = sorted(
                    tuple(
                        set(
                            (
                                germlines[c]["v_gene"][0][0],
                                germlines[c]["j_gene"][0][0],
                            )
                        )
                    )
                )
                germlines_and_antigen[f"AB_{ab_domain}_seq"] = get_sequences(self[c])[c]

            germlines_and_antigen["antigen"] = (
                get_sequences(self.get_antigen()[0])[self.get_antigen()[0].id]
                if len(self.get_antigen()) == 1
                else None
            )
        except Exception as e:
            warnings.warn(
                f"Germline/antigen retrieval failed for {self} with error {str(e)}"
            )

        return germlines_and_antigen

    def get_chain_mapping(self):
        """
        Get a dictionary of chain IDs to chain types.

        Returns:
            dict: Dictionary of chain IDs to chain types
        """
        ab_chain_mapping = {v: k for k, v in self.get_domain_assignment().items()}
        antigen_chain_mapping = {c.id: "Ag" for c in self.get_antigen()}
        chain_mapping = {
            **ab_chain_mapping,
            **antigen_chain_mapping,
        }
        return chain_mapping

    def save(self, save_as=None, ab_only: bool = False, format: str = "pdb"):
        """
        Save AB object as PDB or MMCIF file.

        Args:
            save_as (str, optional): File path to save AB to.
            ab_only (bool, optional): Whether to save antibody only.
            format (str, optional): "pdb" or "mmcif".
        """
        from . import TCRIO

        tcrio = TCRIO.TCRIO()
        tcrio.save(self, save_as=save_as, tcr_only=ab_only, format=format)

    def crop(self, *, remove_het_atoms: bool = True) -> None:
        """
        Crop antibody to variable domain.

        This method mutates the AB object.

        Args:
            remove_het_atoms: remove het atoms from structure as well
        """
        new_child_dict = {}
        for chain in self:
            new_chain = TCRchain(chain.id)  # consider replacing with ABchain

            for residue in chain:
                if residue.id[1] in IMGT_VARIABLE_DOMAIN or (
                    not remove_het_atoms and residue.id[0] != " "
                ):
                    new_chain.add(residue.copy())

            new_chain.analyse(chain.chain_type)
            new_chain.set_engineered(chain.engineered)
            new_chain.xtra.update(chain.xtra)
            new_child_dict[new_chain.id] = new_chain

        for chain_id in new_child_dict:
            del self[chain_id]

        for new_chain in new_child_dict.values():
            self.add(new_chain)

    def standardise_chain_names(self):
        """Raises NotImplementedError."""
        raise NotImplementedError()

    def _validate_chain_standardising(self) -> None:
        if hasattr(self, "antigen") and len(self.antigen) > 1:
            msg = "More than one antigen molecule is not currently supported for standardising."
            raise ABError(msg)

    def _standardise_antigen_chain_names(self) -> None:
        """
        Will give the antigen the chain id C. Does not support multiple antigens.
        """
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonWarning)
            self.antigen[0].id = "C"


class hlAb(Ab):
    """
    Heavy-light antibody class. Inherits from AB.
    This is a subclass of AB for antibodies with heavy and light chains.
    """

    def __init__(self, c1, c2):
        """
        Initialise hlAb object.

        Args:
            c1: heavy or light antibody chain
            c2: heavy or light antibody chain
        """
        if c1.chain_type == "H":
            Entity.__init__(self, c1.id + c2.id)
        else:
            Entity.__init__(self, c2.id + c1.id)

        self.level = "H"
        self._add_domain(c1)
        self._add_domain(c2)
        self.child_list = sorted(
            self.child_list, key=lambda x: x.chain_type, reverse=True
        )  # H -> L
        self.antigen = []
        self.engineered = False
        self.scAb = False  # rare, but possible

    def __repr__(self):
        """
        String representation of the hlAb object.
        """
        return "<AB %s%s heavy=%s; light=%s>" % (self.VH, self.VL, self.VH, self.VL)

    def _add_domain(self, chain):
        """
        Add a variable heavy or variable light domain to the AB object.
        Links the domain to the chain ID.

        Args:
            chain: antibody chain whose domain is being added
        """
        if chain.chain_type == "H":
            self.VH = chain.id
        elif chain.chain_type == "L":
            self.VL = chain.id

        self.add(chain)

    def get_VH(self):
        """
        Retrieve the variable heavy chain of the antibody.

        Returns:
            chain object
        """
        if hasattr(self, "VH"):
            return self.child_dict[self.VH]

    def get_VL(self):
        """
        Retrieve the variable light chain of the antibody.

        Returns:
            chain object
        """
        if hasattr(self, "VL"):
            return self.child_dict[self.VL]

    def get_domain_assignment(self):
        """
        Retrieve the domain assignment of the antibody.

        Returns:
            dict: e.g. {"VH": "H", "VL": "L"}
        """
        try:
            return {"VH": self.VH, "VL": self.VL}
        except AttributeError:
            if hasattr(self, "VH"):
                return {"VH": self.VH}
            if hasattr(self, "VL"):
                return {"VL": self.VL}
        return None

    def is_engineered(self):
        """
        Flag for engineered antibodies.

        Returns:
            bool
        """
        if self.engineered:
            return True
        else:
            vh, vl = self.get_VH(), self.get_VL()
            for var_domain in [vh, vl]:
                if var_domain and var_domain.is_engineered():
                    self.engineered = True
                    return self.engineered

            self.engineered = False
            return False

    def get_fragments(self):
        """
        Retrieve the fragments, i.e. FW and CDR loops of the antibody, as a generator.
        """
        vh, vl = self.get_VH(), self.get_VL()

        for var_domain in [vh, vl]:
            if var_domain:
                for frag in var_domain.get_fragments():
                    yield frag

    def standardise_chain_names(self) -> None:
        """
        Standardise the antibody and antigen chain names to the following convention.

        Convention:
            - C - antigen chain
            - H - antibody heavy chain
            - L - antibody light chain

        Note, this mutates the original object.

        Raises:
            ABError: if there is more than one antigen molecule attached to the antibody.
        """
        self._validate_chain_standardising()

        new_id = []
        new_child_dict = {}

        if hasattr(self, "VH"):
            new_child_dict["H"] = self.child_dict[self.VH]
            self.VH = "H"
            new_id.append("H")

        if hasattr(self, "VL"):
            new_child_dict["L"] = self.child_dict[self.VL]
            self.VL = "L"
            new_id.append("L")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonWarning)

            for chain_id, chain in new_child_dict.items():
                chain.id = chain_id

        self.child_dict = new_child_dict

        if hasattr(self, "antigen") and self.antigen:
            self._standardise_antigen_chain_names()

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonWarning)
            self.id = "".join(new_id)
