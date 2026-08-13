import json
import os
import random
from typing import Any, Dict, List, Optional, Union

import requests

from ..helper.bin2Quaternary import byte2QUATS
from ..helper.quaternary2Bin import quats_to_bytes

# IMPORTANT: if you plan to use this, you should change the MESA_URL to a (local) instance of your own.
MESA_URL = "http://pc12291.mathematik.uni-marburg.de:5000/api/all"


class DNARules2:
    def __init__(self, active_rules=None):
        self.active_rules = []

    @staticmethod
    def sim_mutation(sequences):
        """
        Simulates the mutation of nucleobases with a static probability for a sequence
        :param seq: Sequence to be mutated
        :return: Mutated sequence
        """
        res_seq = []
        for seq in sequences:
            mod_seq = ""
            # randint(0, x) <= y defines the chance to mutate: 1 - y/x * 100% chance to mutate
            for x in seq:
                if random.randint(0, 1000) <= 1000:  # always True
                    mod_seq += x
                else:
                    mod_seq += random.choice(("A", "T", "G", "C"))
            res_seq.append(mod_seq)
        return res_seq

    @staticmethod
    def _decode_mutated_sequence(seq: str) -> bytes:
        dna_data_bin_enc = b""
        for i in range(0, len(seq), 4):
            try:
                dna_data_bin_enc += quats_to_bytes(seq[i : i + 4])
            except ValueError:
                continue
        return dna_data_bin_enc

    @staticmethod
    def _calculate_changes(dna_data_dec: str, seq_org: str) -> int:
        return sum(
            1
            for index in range(min(len(dna_data_dec), len(seq_org)))
            if dna_data_dec[index] != seq_org[index]
        )

    @staticmethod
    def _get_dropchance(seq: str, seq_org: str) -> float:
        dna_data_bin_enc = DNARules2._decode_mutated_sequence(seq)
        dna_data_dec = "".join(byte2QUATS(value) for value in dna_data_bin_enc)
        changes = DNARules2._calculate_changes(dna_data_dec, seq_org)
        if changes >= 3 - abs(len(dna_data_dec) - len(seq_org)):
            return 1.0

        dropchance = 10 * changes
        for index, value in enumerate(dna_data_bin_enc):
            if value != dna_data_bin_enc[index] and dropchance <= 95:
                dropchance += 5
        return dropchance / 100

    @staticmethod
    def apply_all_rules(packet, from_web=True):
        """
        Calculates the dropchance for either a generated packet or a list of generated packets with the option to
        simulate mutations with a static probability or getting the mutated sequences from the website. Also allows
        to import a config file with configurations downloaded from the website.
        :param packet: The packet or list of packets to calculate the dropchance for
        :param from_web: True: Get the mutated sequences from the website. False: Use the static simulation
        #:param from_file: True: Use a config file for the websiterequest. False: Enter configurations manual in @get_mutated
        :return: The dropchance for a packet or a list of packets
        """
        packets = packet if isinstance(packet, list) else [packet]
        dna_data = [pack.get_dna_struct(True) for pack in packets]
        dna_data_mutated = DNARules2.get_mutated(dna_data, from_web)
        res_err = [
            DNARules2._get_dropchance(seq, dna_data[index])
            for index, seq in enumerate(dna_data_mutated)
        ]
        return res_err[0] if len(res_err) == 1 else res_err

    # executes requests to the website with given parameters
    @staticmethod
    def get_mutated_from_web(
        seq: List[str],
        config: Optional[Dict[str, Any]] = None,
        json_config: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        """
        Builds and executes the request for the website from given configurations, either a config file or manual
        added parameters
        :param seq: The sequence(s) to mutate
        :param config: If used, the config file downloaded from the website
        :param json_config: If used, the manually generated json_config
        :return: The response from the website as dictionary
        """
        header = {"content-type": "application/json;charset=UTF-8"}
        payload = dict(
            json_config if json_config is not None else config if config is not None else {}
        )
        payload["sequence"] = seq
        payload["asHTML"] = False
        res = requests.post(MESA_URL, json=payload, headers=header, timeout=30)
        res.raise_for_status()
        response = res.json()
        if not isinstance(response, dict):
            raise TypeError("MESA API response must be a JSON object")
        return response

    @staticmethod
    def get_mutated(
        seq: Union[str, List[str]], from_web: bool, json_config: Optional[Dict[str, Any]] = None
    ) -> List[str]:
        """
        Takes the sequences to mutate and calls either @sim_mutation or @get_mutated_from_web. If the website is used, the
        results are also processed to get only the mutated sequences for the input
        :param seq: The sequence(s) to mutate
        :param from_web: True: Use the website, False: Use the simulated mutation
        :param json_config: If used, the manually generated json_config
        :return: A list of mutated sequences for the input
        """
        seq_list = [seq] if isinstance(seq, str) else seq
        if from_web:
            if json_config is not None:
                res_all = DNARules2.get_mutated_from_web(seq=seq_list, json_config=json_config)
            else:
                try:
                    file = os.environ["dna_sim_config"]
                except KeyError:
                    print("Could not find ENV-Var 'dna_sim_config', falling back to 'mosla.json'")
                    file = "mosla.json"
                with open(file) as json_file:
                    config = json.load(json_file)
                    res_all = DNARules2.get_mutated_from_web(seq=seq_list, config=config)
            res_seq = []
            for sequence in seq_list:
                res_entry = res_all[sequence]
                res_seq.append(res_entry["res"]["modified_sequence"])
            return res_seq
        else:
            res_seq = DNARules2.sim_mutation(seq_list)
        return res_seq


if __name__ == "__main__":
    sequence = "AAAACCCCGGGGTTTT"
    print("Test with the mosla.json file for: " + sequence)
    print(DNARules2.get_mutated(sequence, True))
