"""
Created on 3 April 2024
Adapted from TCRParser to AbParser

AbParser object based on BioPython's PDB parser.
Parses antibody structures, pairs heavy/light chains, and assigns antigens.
"""

from itertools import combinations, product
import sys
import os
from collections import defaultdict
import warnings

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB import NeighborSearch

from .annotate import annotate, extract_sequence, align_numbering
from ..utils.error_stream import ErrorStream

from .AbStructure import AbStructure
from .Model import Model
from .Ab import Ab, hlAb, hkAb
from .Holder import Holder
from .AbChain import AbChain
from .AGchain import AGchain
from Bio.PDB.Residue import Residue
from .Fragment import Fragment
from .Chemical_components import is_aa, is_common_buffer, get_res_type, is_carbohydrate


class AbParser(PDBParser, MMCIFParser):
    def __init__(self, PERMISSIVE=True, get_header=True, QUIET=False):
        """
        Initialise the parser.
        Structures are numbered using ANARCI and IMGT definitions.
        """
        self.pdb_parser = PDBParser(PERMISSIVE, get_header, None, QUIET)
        self.mmcif_parser = MMCIFParser(None, QUIET)
        self.QUIET = QUIET

        self.numbering_method = "anarci"
        self.numbering_scheme = "imgt"
        self.definition = "imgt"

        self.current_file = None

    def _create_chain(self, chain, new_chain_id, numbering, chain_type):
        """
        Create a new antibody chain.
        Residues before the numbered region are ignored.
        """
        if chain_type in ["H", "K", "L"]:
            newchain = AbChain(new_chain_id)
        else:
            raise ValueError(f"Unsupported antibody chain type: {chain_type}")

        newchain.numbering = numbering
        unnumbered_list = []
        added = False

        for residue in chain.get_list():
            add = False
            if residue.id in numbering:
                if numbering[residue.id]:
                    add = True
                    res_id = (
                        residue.id[0],
                        numbering[residue.id][0],
                        numbering[residue.id][1],
                    )

            if add:
                added = True
                newresidue = Residue(res_id, residue.resname, residue.segid)
                for atom in residue.get_list():
                    newresidue.add(atom.copy())
                newresidue.imgt_numbered = True
                newchain.add(newresidue)

            elif added:
                unnumbered_list.append(residue)

        ended = sorted([i for i in numbering.values() if i != ""])[-1][0]
        for residue in unnumbered_list:
            ended += 1
            res_id = (residue.id[0], ended, " ")
            newresidue = Residue(res_id, residue.resname, residue.segid)
            for atom in residue.get_list():
                newresidue.add(atom.copy())
            newchain.add(newresidue)
            newchain.add_unnumbered(newresidue)

        newchain.analyse(chain_type)
        return newchain

    def _create_scAb_chains(
        self, chain, new_chain_id, numbering_1, numbering_2, chain_type1, chain_type2
    ):
        """
        Create two antibody chains from a single physical chain.
        Intended for scFv-like constructs where ANARCI finds two variable domains.
        """
        if chain_type1 != "H" and chain_type2 == "H":
            return self._create_scAb_chains(
                chain,
                new_chain_id,
                numbering_2,
                numbering_1,
                chain_type2,
                chain_type1,
            )

        newchain1 = AbChain(new_chain_id.lower())
        newchain2 = AbChain(new_chain_id.upper())

        newchain1.numbering = numbering_1
        newchain2.numbering = numbering_2
        newchains = [newchain1, newchain2]

        unnumbered_set = set()
        added = False
        numbered_pos = set(numbering_1.keys()) | set(numbering_2.keys())

        for i, numbering in enumerate([numbering_1, numbering_2]):
            newchain = newchains[i]

            for residue in chain.get_list():
                add = False
                if residue.id in numbering:
                    if numbering[residue.id]:
                        add = True
                        res_id = (
                            residue.id[0],
                            numbering[residue.id][0],
                            numbering[residue.id][1],
                        )

                if add:
                    added = True
                    newresidue = Residue(res_id, residue.resname, residue.segid)
                    for atom in residue.get_list():
                        newresidue.add(atom.copy())
                    newresidue.imgt_numbered = True
                    newchain.add(newresidue)

                elif added and residue.id not in numbered_pos:
                    unnumbered_set.add(residue)

        ended = sorted(numbering_1.values())[-1][0]
        for residue in sorted(unnumbered_set, key=lambda z: z.id[1]):
            ended += 1
            res_id = (residue.id[0], ended, " ")
            newresidue = Residue(res_id, residue.resname, residue.segid)
            for atom in residue.get_list():
                newresidue.add(atom.copy())

            newchain1.add(newresidue)
            newchain1.add_unnumbered(newresidue)

        newchain1.analyse(chain_type1)
        newchain2.analyse(chain_type2)

        return newchain1, newchain2

    def _number_and_annotate_chain(self, chain, prenumbering=None, ali_dict={}):
        germline_info = False

        if prenumbering and chain.id in prenumbering:
            if len(prenumbering[chain.id]) == 2:
                numbering = [{}, {}]
                region_types = ["", ""]

                numbering[0], region_types[0] = self._prenumbered(
                    chain, prenumbering, ali_dict, n=0
                )
                numbering[1], region_types[1] = self._prenumbered(
                    chain, prenumbering, ali_dict, n=1
                )

                rtypes = set(region_types)

                if rtypes == {"H", "K"} or rtypes == {"H", "L"}:
                    chain_type = "".join(region_types)
                    scAb = True
                else:
                    chain_type = region_types[0]
                    numbering = numbering[0]
                    scAb = False
                    print(
                        "Warning multiple variable regions of unsupported type(s) "
                        f"({region_types}) found on chain {chain.id}. "
                        "Taking the first variable region only.",
                        file=self.warnings,
                    )

            elif prenumbering[chain.id][0][-1] not in ["H", "K", "L"]:
                numbering, chain_type, germline_info, scAb = annotate(chain)

            else:
                numbering, chain_type = self._prenumbered(
                    chain, prenumbering, ali_dict, n=0
                )
                scAb = False

        else:
            numbering, chain_type, germline_info, scAb = annotate(chain)

        return numbering, chain_type, germline_info, scAb

    def _get_header_info(self, abstructure, chain, germline_info):
        if chain.id in abstructure.header["chain_details"]:
            engineered = abstructure.header["chain_details"][chain.id]["engineered"]
            details = abstructure.header["chain_details"][chain.id]
        else:
            engineered = False
            details = {"molecule": "unknown", "engineered": False}

        details["genetic_origin"] = germline_info
        return details, engineered

    def _read_structure_file(self, file, id):
        _, ext = os.path.splitext(file)
        if ext.lower() == ".pdb":
            structure = self.pdb_parser.get_structure(id, file)
            self.current_parser = self.pdb_parser
        elif ext.lower() in [".cif", ".mmcif"]:
            structure = self.mmcif_parser.get_structure(id, file)
            self.current_parser = self.mmcif_parser
        else:
            self.warnings.write(f"Unrecognised structure file format: {file}")
            raise ValueError

        abstructure = AbStructure(structure.id)
        abstructure.set_header(structure.header)
        self._analyse_header(abstructure)
        return structure, abstructure

    def _initialise_model(self, model):
        newmodel = Model(model.id)

        agchains = Holder("Antigen")
        abchains = Holder("AbChain")

        newmodel.add(agchains)
        newmodel.add(abchains)
        return newmodel, agchains, abchains

    def _instantiate_ab(self, chain1, chain2):
        obs_chaintypes = {chain1.chain_type, chain2.chain_type}

        if obs_chaintypes == {"H", "L"}:
            return hlAb(chain1, chain2)
        elif obs_chaintypes == {"H", "K"}:
            return hkAb(chain1, chain2)
        return None

    def get_ab_structure(
        self,
        id,
        file,
        prenumbering=None,
        ali_dict={},
        crystal_contacts=[],
        include_symmetry_mates=True,
    ):
        """
        Post-processing of a Bio.PDB structure object into an antibody context.

        Args:
            id: identifier for the structure
            file: path to the input structure file

        Optional:
            prenumbering: prenumbering information for chains in the structure
        """
        self.warnings = ErrorStream()
        self.include_symmetry_mates = include_symmetry_mates
        self.current_file = file

        structure, abstructure = self._read_structure_file(file, id)

        for mid in range(len(structure.child_list) - 1, -1, -1):
            model = structure.child_list[mid]
            newmodel, agchains, abchains = self._initialise_model(model)
            abstructure.add(newmodel)

            for chain in model.get_list():
                numbering, chain_type, germline_info, scAb = (
                    self._number_and_annotate_chain(chain, prenumbering, ali_dict)
                )

                details, engineered = self._get_header_info(
                    abstructure, chain, germline_info
                )

                if numbering and chain_type in ["H", "K", "L"]:
                    newchain = self._create_chain(
                        chain, chain.id, numbering, chain_type
                    )
                    newchain.set_engineered(engineered)
                    newchain.xtra.update(details)
                    abchains.add(newchain)

                elif numbering and scAb:
                    types = list(chain_type)
                    domain1, domain2 = numbering

                    chain1, chain2 = self._create_scAb_chains(
                        chain, chain.id, domain1, domain2, types[0], types[1]
                    )
                    chain1.set_engineered(engineered)
                    chain1.xtra.update(details)
                    chain2.set_engineered(engineered)
                    chain2.xtra.update(details)

                    try:
                        if (
                            chain1.child_dict[(" ", 104, " ")]["CA"]
                            - chain2.child_dict[(" ", 104, " ")]["CA"]
                        ) <= 22:
                            ab = self._instantiate_ab(chain1, chain2)
                            if ab is not None:
                                ab.scAb = True
                                newmodel.add(ab)
                                if chain1.id in abchains:
                                    abchains.detach_child(chain1.id)
                                if chain2.id in abchains:
                                    abchains.detach_child(chain2.id)
                            else:
                                abchains.add(chain1)
                                abchains.add(chain2)
                        else:
                            abchains.add(chain1)
                            abchains.add(chain2)
                    except KeyError:
                        abchains.add(chain1)
                        abchains.add(chain2)

                else:
                    newchain = self._create_ag_chain(chain)
                    newchain.set_engineered(engineered)
                    newchain.xtra.update(details)
                    agchains.add(newchain)

            pairings = self._pair_chains(abchains)
            for pair in pairings:
                if pair[0].id in abchains:
                    abchains.detach_child(pair[0].id)
                if pair[1].id in abchains:
                    abchains.detach_child(pair[1].id)

                ab = self._instantiate_ab(pair[0], pair[1])
                if ab is None:
                    self.warnings.write(
                        "Unusual pairing between %s (V%s) and %s (V%s) has been detected. "
                        "Treating as separate antibody chains.\n"
                        % (
                            pair[0].id,
                            pair[0].chain_type,
                            pair[1].id,
                            pair[1].chain_type,
                        )
                    )
                    abchains.add(pair[0])
                    abchains.add(pair[1])
                    continue

                newmodel.add(ab)

            self._match_units(newmodel, abchains, agchains, crystal_contacts)
            del structure.child_list[mid]

            empty_holders = [
                holder.id for holder in newmodel.child_list if not holder.child_list
            ]
            for holder_id in empty_holders:
                newmodel.detach_child(holder_id)

        del structure
        if not self.QUIET and self.warnings.log:
            sys.stderr.write("\n".join(self.warnings.log))
            sys.stderr.write("\n")
        abstructure.warnings = self.warnings

        self.current_file = None
        return abstructure

    def _analyse_header(self, header):
        """
        Analysis of the header parsed by Biopython.
        Adds per-chain metadata such as molecule and engineered flags.
        """
        if isinstance(header, AbStructure):
            header = header.get_header()
        elif not header:
            header = {}

        header["chain_details"] = {}
        if "compound" in header:
            for compound in header["compound"]:
                if "chain" in header["compound"][compound]:
                    chains = [
                        c.strip().upper()
                        for c in header["compound"][compound]["chain"].split(",")
                        if len(c.strip()) == 1
                    ]

                    for chain in chains:
                        if chain not in header["chain_details"]:
                            header["chain_details"][chain] = {}

                    if "molecule" in header["compound"][compound]:
                        for chain in chains:
                            header["chain_details"][chain]["molecule"] = header[
                                "compound"
                            ][compound]["molecule"]
                    else:
                        for chain in chains:
                            header["chain_details"][chain]["molecule"] = "unknown"

                    if "engineered" in header["compound"][compound]:
                        if (
                            "no" in header["compound"][compound]["engineered"]
                            or "false" in header["compound"][compound]["engineered"]
                            or not header["compound"][compound]["engineered"]
                        ):
                            header["compound"][compound]["engineered"] = False
                        else:
                            header["compound"][compound]["engineered"] = True
                        for chain in chains:
                            header["chain_details"][chain]["engineered"] = header[
                                "compound"
                            ][compound]["engineered"]
                    else:
                        for chain in chains:
                            header["chain_details"][chain]["engineered"] = False
                else:
                    continue
        else:
            sys.stderr.write("Header could not be parsed")

    def _create_ag_chain(self, chain):
        """
        Create a new antigen chain.
        This simply means a chain that is not annotated as an antibody chain.
        """
        newchain = AGchain(chain.id)
        for residue in chain.get_list():
            newresidue = Residue(residue.id, residue.resname, residue.segid)
            newchain.add(newresidue)
            for atom in residue.get_list():
                newresidue.add(atom.copy())
        newchain.set_type()
        return newchain

    def _pair_chains(self, chains):
        """
        Pair heavy-light and heavy-kappa chains to form antibodies.
        Uses the CA atom at IMGT 104 as a simple interface heuristic.
        Greedily selects the closest non-overlapping valid pairs.
        """
        candidates = []
        points = {
            "H": (" ", 104, " "),
            "K": (" ", 104, " "),
            "L": (" ", 104, " "),
        }

        for pair in combinations(chains, 2):
            types = {pair[0].chain_type, pair[1].chain_type}
            if types not in [{"H", "K"}, {"H", "L"}]:
                continue

            try:
                a1 = pair[0].child_dict[points[pair[0].chain_type]]["CA"]
                a2 = pair[1].child_dict[points[pair[1].chain_type]]["CA"]
            except KeyError:
                continue

            dist = a1 - a2
            if dist < 22:
                candidates.append((dist, pair))

        candidates.sort(key=lambda x: x[0])

        used = set()
        pairings = []
        for _, pair in candidates:
            if pair[0].id in used or pair[1].id in used:
                continue
            used.add(pair[0].id)
            used.add(pair[1].id)
            pairings.append(pair)

        return pairings

    def _get_sugar_fragments(self, sugar):
        """
        Get connected hetatoms to form sugar molecules.
        """
        sugar = dict(list(zip([s.id for s in sugar], sugar)))

        connect_records = {}
        for c in [
            line.strip() for line in self.current_parser.trailer if "CONECT" in line
        ]:
            try:
                connect_records[int(c[6:11])] = []
            except IndexError:
                continue
            for b, e in [(11, 16), (16, 21), (21, 26), (26, 31)]:
                try:
                    if c[b:e].strip():
                        connect_records[int(c[6:11])].append(int(c[b:e]))
                    else:
                        break
                except IndexError:
                    break
                except ValueError:
                    self.warnings.write(
                        "Warning: unexpected CONECT record format %s" % c.strip()
                    )

        monomer_atoms = []
        polymers = []
        if connect_records:
            atomid_to_resid = {}
            for r in sugar:
                for atom in sugar[r]:
                    atomid_to_resid[atom.serial_number] = sugar[r].id

            r_connections = {}
            for a in connect_records:
                if a in atomid_to_resid:
                    try:
                        r_connections[atomid_to_resid[a]].update(
                            [
                                atomid_to_resid[ai]
                                for ai in connect_records[a]
                                if ai in atomid_to_resid
                            ]
                        )
                    except KeyError:
                        r_connections[atomid_to_resid[a]] = set(
                            [
                                atomid_to_resid[ai]
                                for ai in connect_records[a]
                                if ai in atomid_to_resid
                            ]
                        )

            connected_sets = []
            for r in sorted(r_connections, key=lambda x: x[1]):
                added = 0
                for i in range(len(connected_sets)):
                    if connected_sets[i] & r_connections[r]:
                        connected_sets[i].update(r_connections[r])
                        added = 1
                        break
                if not added:
                    connected_sets.append(r_connections[r])

            n = 0
            for mol in connected_sets:
                if len(mol) > 1:
                    polymers.append(Fragment("sugar%d" % n))
                    for r in sorted(mol, key=lambda x: x[1]):
                        polymers[n].add(sugar[r])
                    n += 1
                else:
                    monomer_atoms += [atom for atom in sugar[list(mol)[0]]]

        else:
            for s in sugar:
                monomer_atoms += [atom for atom in sugar[s]]

        return polymers, monomer_atoms

    def _find_chain_hetatoms(self, chain):
        """
        Filter out HETATM records from antibody chains for downstream antigen assignment.
        """
        hetatoms, sugars = [], []
        for residue in chain.get_unnumbered():
            if residue.id[0] == "W" or is_aa(residue, standard=False):
                continue
            if is_carbohydrate(residue):
                sugars.append(residue)
            else:
                hetatoms.extend(list(residue.get_atoms()))

        return hetatoms, sugars

    def _prepare_ab(self, ab, cdr_atoms, antigen_hetatoms, antigen_sugars):
        """
        Prepare antibody contact atoms and chain-associated HETATMs.
        For antibodies, all CDRs are used, not just CDR3.
        """
        for cdr in ab.get_CDRs():
            cdr_atoms[ab.id] += [
                atom for atom in cdr.get_atoms() if atom.id == "CB" or atom.id == "CA"
            ]

        if isinstance(ab, Ab):
            for chain in ab.get_chains():
                antigen_hetatoms[chain.id], antigen_sugars[chain.id] = (
                    self._find_chain_hetatoms(chain)
                )

        elif isinstance(ab, AbChain):
            antigen_hetatoms[ab.id], antigen_sugars[ab.id] = self._find_chain_hetatoms(
                ab
            )

    def _prepare_abs_and_antigens_for_pairing(
        self,
        model,
        antibodies,
        agchains,
        crystal_contacts,
    ):
        antigen_atoms, cdr_atoms, antigen_hetatoms, antigen_sugars = (
            defaultdict(list),
            defaultdict(list),
            defaultdict(list),
            defaultdict(list),
        )

        for ab in antibodies:
            self._prepare_ab(ab, cdr_atoms, antigen_hetatoms, antigen_sugars)

        for antigen in agchains:
            antigen_atoms[antigen.id] = [
                a
                for a in antigen.get_atoms()
                if a.parent.id[0] == " " or is_aa(a.parent)
            ]
            antigen_hetatoms[antigen.id] = [
                a
                for a in antigen.get_atoms()
                if a.parent.id[0].startswith("H") and not is_aa(a.parent)
            ]

        sugars = []
        for chain_id in antigen_sugars:
            if antigen_sugars[chain_id]:
                polymers, monomer_atoms = self._get_sugar_fragments(
                    antigen_sugars[chain_id]
                )
                sugars += polymers
                antigen_hetatoms[chain_id] += monomer_atoms

        non_empty_ag = [k for k in antigen_hetatoms if antigen_hetatoms[k]]

        self._protein_peptide_pass(
            model, antibodies, cdr_atoms, antigen_atoms, crystal_contacts
        )
        self._het_sugar_pass(
            antibodies,
            cdr_atoms,
            non_empty_ag,
            antigen_hetatoms,
            sugars,
            distance=8.0,
        )

        return (
            model,
            antibodies,
            agchains,
            crystal_contacts,
            antigen_atoms,
            cdr_atoms,
            antigen_hetatoms,
            antigen_sugars,
        )

    def _has_polymeric_antigen_candidates(self, agchains):
        """
        Only trigger symmetry-mate search if another polymeric chain exists
        that could plausibly be an antigen.
        """
        for chain in agchains:
            if hasattr(chain, "get_type") and chain.get_type() in {
                "protein",
                "peptide",
            }:
                return True
        return False

    def _add_symmetry_mate_antigens(self, symmetry_mates, agchains):
        """
        Add antigen-holder chains from symmetry mates into the current model's
        antigen holder as extra candidate antigens.
        """
        for mate in symmetry_mates or []:
            for holder in mate.get_holders():
                if holder.id != "Antigen":
                    continue
                for ag in holder:
                    if ag.id not in agchains.child_dict:
                        agchains.add(ag.copy())

    def _match_units(self, model, abchains, agchains, crystal_contacts=[]):
        """
        Match antigen chains / ligands to antibodies.
        """
        antibodies = [h for h in model if isinstance(h, Ab)] + abchains.child_list

        (
            model,
            antibodies,
            agchains,
            crystal_contacts,
            antigen_atoms,
            cdr_atoms,
            antigen_hetatoms,
            antigen_sugars,
        ) = self._prepare_abs_and_antigens_for_pairing(
            model,
            antibodies,
            agchains,
            crystal_contacts,
        )

        unpaired_antibody_exists = any(
            not getattr(ab, "antigen", []) for ab in antibodies
        )

        if (
            self.include_symmetry_mates
            and unpaired_antibody_exists
            and self._has_polymeric_antigen_candidates(agchains)
        ):
            try:
                symmetry_mates = self._generate_symmetry_mates()
                self._add_symmetry_mate_antigens(symmetry_mates, agchains)
            except Exception as e:
                warnings.warn(f"Symmetry mate generation failed with: {str(e)}")

            (
                model,
                antibodies,
                agchains,
                crystal_contacts,
                antigen_atoms,
                cdr_atoms,
                antigen_hetatoms,
                antigen_sugars,
            ) = self._prepare_abs_and_antigens_for_pairing(
                model,
                antibodies,
                agchains,
                crystal_contacts,
            )

    def _generate_symmetry_mates(self):
        print("Generating symmetry mates to pair antigens.")
        from .utils.symmetry_mates import (
            get_symmetry_mates,
        )

        return get_symmetry_mates(self.current_file)

    def _protein_peptide_pass(
        self, model, complexes, receptor_atoms, antigen_atoms, crystal_contacts=[]
    ):
        """
        Generic method to assign protein/peptide antigens to antibody objects.
        """
        ns = NeighborSearch(
            [atom for chain in receptor_atoms for atom in receptor_atoms[chain]]
            + [atom for chain in antigen_atoms for atom in antigen_atoms[chain]]
        )
        contacts = [con for con in ns.search_all(8.0, "R")]
        contact_freq = defaultdict(lambda: defaultdict(int))

        all_cpx_chains = {}
        for cpx in complexes:
            cpx_ch = list(cpx.id)
            for c in cpx_ch:
                all_cpx_chains[c] = cpx.id

        cpxids = set(all_cpx_chains.values())
        ags = set()

        for c in contacts:
            p1 = str(c[0].parent.id)
            p2 = str(c[1].parent.id)

            potential_contact = p1 + p2
            potential_contact2 = p2 + p1

            if (
                p1 == p2
                or potential_contact in contact_freq
                or potential_contact2 in contact_freq
            ):
                continue

            if (
                (potential_contact not in cpxids)
                and (p1 in all_cpx_chains)
                and (p2 not in all_cpx_chains)
            ):
                C = all_cpx_chains[p1]
                ag = p2
            elif (
                (potential_contact2 not in cpxids)
                and (p2 in all_cpx_chains)
                and (p1 not in all_cpx_chains)
            ):
                C = all_cpx_chains[p2]
                ag = p1
            else:
                continue

            contact_freq[C][ag] += 1
            ags.add(ag)

        for cpx_id in cpxids:
            if contact_freq[cpx_id]:
                ag = max(contact_freq[cpx_id], key=lambda x: contact_freq[cpx_id][x])

                if (cpx_id, ag) not in crystal_contacts:
                    model[cpx_id].antigen = []
                    model[cpx_id]._add_antigen(model[ag])

                    if ag in ags:
                        ags.remove(ag)

        for ag in ags:
            cmax = 0
            for C in contact_freq:
                if ag in contact_freq[C] and (C, ag) not in crystal_contacts:
                    if contact_freq[C][ag] > cmax:
                        paired_cpx = C
                        cmax = contact_freq[C][ag]
            if cmax:
                if len(contact_freq) > 1:
                    self.warnings.write(
                        "Crystal Contact Warning: antigen %s has been paired with antibody %s"
                        % (str(ag), str(paired_cpx))
                    )
                    model[paired_cpx]._add_antigen(model[ag])
                else:
                    model[paired_cpx]._add_antigen(model[ag])

    def _het_sugar_pass(
        self,
        receptors,
        receptor_atoms,
        non_empty_ag,
        antigen_hetatoms,
        sugars,
        distance=8.0,
    ):
        """
        Assign nearby HETATM/sugar antigens to antibodies.
        """
        for rec, antigen_het in product(receptors, non_empty_ag):
            ns = NeighborSearch(antigen_hetatoms[antigen_het])

            for atom in receptor_atoms[rec.id]:
                contacts = ns.search(atom.get_coord(), distance, level="R")
                if contacts:
                    for contact in contacts:
                        if self._check_het_antigen(contact):
                            residue_type = get_res_type(contact)

                            if residue_type == "Hapten":
                                self.warnings.write(
                                    """Warning: Multiple hapten-antigen like molecules found in binding site -
                                    this needs attention as could be solvent/cofactor."""
                                )
                            if residue_type == "non-polymer":
                                contact.type = "Hapten"
                                contact.get_type = lambda: "Hapten"
                            elif residue_type == "nucleic-acid":
                                contact.type = "nucleic-acid"
                                contact.get_type = lambda: "nucleic-acid"
                            elif residue_type == "saccharide":
                                contact.type = "carbohydrate"
                                contact.get_type = lambda: "carbohydrate"
                            rec._add_antigen(contact)

    def _check_het_antigen(self, residue):
        """
        Perform checks on a potential HETATM antigen residue.
        """
        if is_aa(residue, standard=False):
            return False

        if is_common_buffer(residue):
            if not self.QUIET:
                self.warnings.write(
                    "Common molecule %s found in the binding site - not considered an antigen"
                    % residue.get_resname()
                )
            return False

        return True

    def _prenumbered(self, chain, prenumbering, ali_dict={}, n=0):
        """
        Deal with numbering supplied by the user.
        """
        if ali_dict:
            ali_dict = ali_dict[chain.id][n]

        annotation, chain_type = prenumbering[chain.id][n]

        try:
            sequence_list, sequence_str, warnings = extract_sequence(
                chain, return_warnings=True
            )
            numbering = align_numbering(annotation, sequence_list, ali_dict)
        except AssertionError:
            sequence_list, sequence_str, warnings = extract_sequence(
                chain, return_warnings=True, ignore_hets=True
            )
            numbering = align_numbering(annotation, sequence_list, ali_dict)

        self.warnings.log += warnings
        return numbering, chain_type
