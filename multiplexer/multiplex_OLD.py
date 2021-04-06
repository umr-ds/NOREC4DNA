import multiprocessing
import random
from functools import partial
from math import ceil, floor
from norec4dna.helper.RU10Helper import intermediate_symbols, from_true_false_list
from norec4dna import RU10Encoder, Encoder, nocode, RU10Decoder, reed_solomon_encode, crc32
from norec4dna.distributions.RaptorDistribution import RaptorDistribution
from norec4dna.rules.FastDNARules import FastDNARules
import numpy as np
from Cryptodome.Protocol.SecretSharing import Shamir
import struct
import secrets
from bitarray import bitarray
from norec4dna.Packet import ParallelPacket
from helpful_scripts.automatedfindminimum import AutomatedFindMinimum


class Multiplexer:
    def __init__(self, file, no_channel, factor=0.75, pack_factor=None, chunk_size=50, as_dna=True, fast=False,
                 secure=False, no_excl=1, id_len_format='H'):
        self.file = file
        self.no_channel = no_channel
        self.factor = factor
        self.pack_factor = factor
        self.id_len_format = id_len_format
        self.chunk_size = chunk_size
        self.as_dna = as_dna
        self.encoder, self.decoder = self.init_multiplex()
        self.no_chunks = self.encoder.number_of_chunks
        self.channel = self.init_channel()
        self.fast = fast
        self.secure = secure
        self.no_excl = no_excl
        if self.fast and self.secure:
            print("Please choose only one option from fast and secure.")
        # prepare the sets of chunks per channel
        if not self.fast and not self.secure:
            self.alld_chunks, self.excl_chunks = self.prep_chunk_sets(self.no_chunks, self.no_channel, self.factor)
        elif self.fast:
            self.alld_chunks, self.excl_chunks = self.prep_chunk_sets_fast(self.no_chunks, self.no_channel,
                                                                           self.no_excl)
        elif self.secure:
            if self.no_channel < 2:
                print("Secure mode is only possible for 2 or more channel.")
            self.used_chunks = [[] for _ in range(0, no_channel)]
            self.alld_chunks, self.excl_chunks = self.prep_chunk_sets_fast(self.no_chunks, self.no_channel,
                                                                           self.no_excl)

    def init_multiplex(self):
        """
        Prepares the packet generating encoder and the overall decoder which needs to be able to decode.
        :return:
        """
        no_chunks = Encoder.get_number_of_chunks_for_file_with_chunk_size(self.file, self.chunk_size,
                                                                          insert_header=False)
        dist = RaptorDistribution(no_chunks)
        if self.as_dna:
            rules = FastDNARules()
        else:
            rules = None
        enc = RU10Encoder(self.file, no_chunks, dist, chunk_size=self.chunk_size, insert_header=False, rules=rules,
                          error_correction=nocode, id_len_format=self.id_len_format, number_of_chunks_len_format="B",
                          save_number_of_chunks_in_packet=False, mode_1_bmp=False)
        enc.prepare()
        dec = RU10Decoder.pseudo_decoder(enc.number_of_chunks, False)
        if dec.distribution is None:
            dec.distribution = RaptorDistribution(enc.number_of_chunks)
            dec.number_of_chunks = enc.number_of_chunks
            _, dec.s, dec.h = intermediate_symbols(enc.number_of_chunks, dec.distribution)
            dec.createAuxBlocks()
        return enc, dec

    def init_channel(self):
        """
        Prepares the channel to send the packets to. Atm these are pseudodecoder. Basically it's a list of objects that
        have an inputNewPacket() method implemented.
        :return:
        """
        channel = [RU10Decoder.pseudo_decoder(self.no_chunks, False) for _ in range(0, self.no_channel)]
        return channel

    @staticmethod
    def prep_chunk_sets_fast(no_chunks, no_channel, no_excl):
        """

        :param no_chunks:
        :param no_channel:
        :param no_excl:
        :return:
        """
        alld_sets = [set() for _ in range(0, no_channel)]
        excl_sets = [set() for _ in range(0, no_channel)]
        chunk_lst = [x for x in range(0, no_chunks)]
        for x in range(0, no_channel):
            excl_sets[x].update(np.random.choice(chunk_lst, no_excl, replace=False))
            for ch in excl_sets[x]:
                alld_sets[x].add(ch)
                chunk_lst.remove(ch)
        for x in range(0, no_channel):
            alld_sets[x].update(chunk_lst)
        return alld_sets, excl_sets

    @staticmethod
    def prep_chunk_sets(no_chunks, no_channel, fac):
        """
        Prepares sets of allowed and exclusive chunks for the given number of channels. The chunks will be distributed
        evenly among the channels. Afterwards a part (1-fac) of these chunks per channel will be declared as exclusive,
        which means they won't be send through another channel. To make decoding possible intersections will be generated.
        Every channels gets some chunks that are allowed in other channels until every channel has 1+fac * chunks/channels
        allowed and 1-fac * chunks/channels exclusive chunks.
        :param no_chunks:
        :param no_channel:
        :param fac:
        :return:
        """
        alld_sets = [set() for _ in range(0, no_channel)]
        excl_sets = [set() for _ in range(0, no_channel)]
        no_add_inters = no_chunks / no_channel * fac
        no_excl = round(no_chunks / no_channel * (1 - fac))
        chunk_lst = [x for x in range(0, no_chunks + 1)]
        chunk_set = set(chunk_lst)
        for x in range(0, no_channel):
            if x == no_channel - 1:  # fill last set
                alld_sets[x].update(chunk_lst)
            else:
                alld_sets[x].update(np.random.choice(chunk_lst, round(no_chunks / no_channel), replace=False))
            for ch in alld_sets[x]:
                chunk_lst.remove(ch)
            excl_sets[x].update(random.sample(alld_sets[x], no_excl))  # set exclusive chunks
        all_excl = set()
        for ex in excl_sets:
            all_excl |= ex
        allwd_chunks = chunk_set - all_excl
        for x in range(0, no_channel):
            alld_sets[x].update(random.sample(allwd_chunks - alld_sets[x], round(no_add_inters)))
        return alld_sets, excl_sets

    def do_multiplex(self, pck_fac=None):
        """
        Fac determines the intersection and exclusive chunks alike.
        1 - Fac is the proportion of chunks that are exclusive for every channel.
        1 + Fac is the proportion of chunks that are allowed for every channel.
        Pck_fac determines the number of packets to generate. If none packets will be generated until decoding is possible.
        Number of packets to generate = no_channel * no_chunks * pck_fac
        :param pck_fac:
        :return:
        """
        used_pcks = [[] for _ in self.channel]
        used_all = 0
        if pck_fac is not None:  # fix amount of packets
            pack_cnt = int(self.no_channel * self.no_chunks * pck_fac)
            for _ in range(0, pack_cnt):
                res = self.create_sort_pack()
                if res is not None:
                    used_pcks[res[0]].append(res[1])
                    used_all += 1
            print(
                "Decoded: " + str(self.decoder.is_decoded()) + " (with fixed packet count). Generated packets: " + str(
                    pack_cnt) + ". Used packets: " + str(used_all) + ". Used " + str(
                    round(used_all / pack_cnt * 100)) + "%.")
        else:  # packets until decoding is possible
            gen_pcks = 0
            while self.decoder.GEPP is None or not self.decoder.is_decoded():
                res = self.create_sort_pack()
                if res is not None:
                    used_pcks[res[0]].append(res[1])
                    used_all += 1
                gen_pcks += 1
            print("Decoded: " + str(self.decoder.is_decoded()) + ". Generated packets: " + str(
                gen_pcks) + ". Used packets: " + str(used_all) + ". Used " + str(
                round(used_all / gen_pcks * 100)) + "%.")
        return used_pcks

    def do_multiplex_with_packets(self, packet_lst, pseudo_decode=False):
        """
        Takes a list of pre generated packets and distributes them to the channels of the multiplexer with the option to
        use pseudo decoder to make sure no channel is able to decode the file encoded in the packets.
        :param packet_lst:
        :param pseudo_decode:
        :return:
        """
        used_packs = [0 for _ in self.channel]
        distribution = RaptorDistribution(self.no_chunks)
        for packet in packet_lst:
            if type(packet) is ParallelPacket:
                packet = AutomatedFindMinimum.parallel_to_normal(packet, nocode, distribution)
            pack_used_chunks = from_true_false_list(self.decoder.removeAndXorAuxPackets(packet))
            if self.secure:
                ch = self.choose_channel_sec(self.alld_chunks, pack_used_chunks, self.used_chunks)
            else:
                ch = self.choose_channel(self.alld_chunks, pack_used_chunks)
            if ch is not None:
                if self.secure:
                    self.used_chunks[ch].append(pack_used_chunks)
                if pseudo_decode:
                    self.decoder.inputNewPacket(packet)
                    self.channel[ch].inputNewPacket(packet)
                used_packs[ch] += 1
        if pseudo_decode:
            print("All packets are distributed. Trying to decode the file.")
        for x in range(0, self.no_channel):
            print("Channel " + str(x) + " used packets: " + str(used_packs[x]) + ".")
            if pseudo_decode and self.channel[x].GEPP is not None and self.channel[x].is_decoded():
                print("Channel " + str(x) + " was able to decode. Simulation failed!")
            else:
                print("Channel " + str(x) + " wasn't able to decode.")
        print("Overall used packets: " + str(sum(used_packs)) + "/" + str(len(packet_lst)) + ": " + str(
            round(sum(used_packs) / len(packet_lst), 4) * 100) + "%.")
        if pseudo_decode and self.decoder.GEPP is not None or self.decoder.is_decoded():
            print("The overall decoder was able to decode the file.")

    def create_sort_pack(self, input_pck=True):
        """
        Creates a new packet, calls choose_channel(_sec) and puts the packet in the right channel. Also returns the
        channel number and the packet.
        :param input_pck:
        :return:
        """
        pack = self.encoder.create_new_packet()
        pack_used_chunks = from_true_false_list(self.decoder.removeAndXorAuxPackets(pack))
        if self.secure:
            ch = self.choose_channel_sec(self.alld_chunks, pack_used_chunks, self.used_chunks)
        else:
            ch = self.choose_channel(self.alld_chunks, pack_used_chunks)  # choose channel for the packet
        if ch is not None:
            if input_pck:
                self.decoder.inputNewPacket(pack)
                self.channel[ch].inputNewPacket(pack)
                if self.secure:
                    self.used_chunks[ch].append(pack_used_chunks)
                if self.channel[ch].GEPP is not None and self.channel[ch].is_decoded():
                    print("Failed. At least one channel was able to decode.")
            return ch, pack
        return None

    @staticmethod
    def choose_channel(alld_chunks, pack_used_chunks):
        """
        Chooses channel for packet based on the allowed chunks of the channel and the chunks in the packet.
        :param alld_chunks:
        :param pack_used_chunks:
        :return:
        """
        for ind, ch in enumerate(alld_chunks):
            if all([x in ch for x in pack_used_chunks]):
                return ind
        return None

    @staticmethod
    def choose_channel_sec(alld_chunks, pack_used_chunks, used_chunks):
        """
        Chooses channel for packet based on the allowed chunks of the channel and the chunks in the packet. Also only
        sends packets with even/uneven grades to channels with even/uneven indices.
        :param alld_chunks:
        :param pack_used_chunks:
        :param used_chunks:
        :return:
        """
        if len(pack_used_chunks) > 1:
            for ind, ch in enumerate(alld_chunks):
                if (ind % 2 == len(pack_used_chunks) % 2 == 0) and all([x in ch for x in pack_used_chunks]):
                    return ind
                elif (ind % 2 == len(pack_used_chunks) % 2 == 1) and all([x in ch for x in pack_used_chunks]):
                    if len(used_chunks[ind]) == 0:
                        return ind
                    for s in used_chunks[ind]:
                        if not set(pack_used_chunks).isdisjoint(set(s)):
                            return None
                    return ind
        return None

    def pack_for_channel(self, channel, no_packs, input_pck=False):
        """
        Creates no_packs packets for the channel with the given number.
        :param channel: 0 - no_channel
        :param no_packs:
        :param input_pck: Put packets into pseudo decoder
        :return:
        """
        pcks = []
        while len(pcks) <= no_packs:
            res = self.create_sort_pack(input_pck=input_pck)
            if res is not None and res[0] == channel:
                pcks.append(res[1])
        return pcks

    def get_noc_shares(self, min):
        """
        Get shares based on the number of active channels.
        :param min:
        :return:
        """
        byte_str = self.no_chunks.to_bytes(16, byteorder="little")
        act_chans = sum([1 for x in self.alld_chunks if len(x) > 0])
        if act_chans < min:
            print("Please decrease the minimum since there are less channels.")
            return None
        return Shamir.split(min, act_chans, byte_str)

    def get_noc_shares_with_data(self, min):
        """
        16 Byte secret:
        number of chunks 4 byte + use_header and err_correction encoding 3 bit + 5 0 bits = 1 byte + 9 random bytes
        :param min:
        :return:
        """
        act_chans = sum([1 for x in self.alld_chunks if len(x) > 0])
        if act_chans < min:
            print("Please decrease the minimum since there are less channels.")
            return None
        noc = struct.pack('I', self.no_chunks)
        arr = bitarray()
        arr.append(self.encoder.insert_header)
        if self.encoder.error_correction is nocode:  # 00 = nocode
            arr.extend([False, False])
        elif self.encoder.error_correction is crc32:  # 01 = crc
            arr.extend([False, True])
        elif self.encoder.error_correction is reed_solomon_encode:  # 10 = reedsolomon
            arr.extend([True, False])
        while (len(arr) < 8):
            arr.append(random.getrandbits(1))
        fill = secrets.token_bytes(11)
        secret = noc + bytes(arr) + fill
        return Shamir.split(min, act_chans, secret)

    @staticmethod
    def combine_noc_shares(shares):
        """
        Combine shares to get the number of chunks back.
        :param shares:
        :return:
        """
        byte_str = Shamir.combine(shares)

        return int.from_bytes(byte_str, byteorder="little")

    @staticmethod
    def combine_noc_shares_with_data(shares):
        byte_str = Shamir.combine(shares)
        noc = struct.unpack('I', byte_str[:4])[0]
        bool_byte = byte_str[4:5]
        insert_header = bool((bool_byte[0] >> 7) & 1)
        err_cor_1 = bool((bool_byte[0] >> 6) & 1)
        err_cor_2 = bool((bool_byte[0] >> 5) & 1)
        if not err_cor_1 and not err_cor_2:
            error_correction = nocode
        elif err_cor_1 is False and err_cor_2 is True:
            error_correction = crc32
        elif err_cor_1 is True and err_cor_2 is False:
            error_correction = reed_solomon_encode
        else:
            raise RuntimeError("This should never happen:", err_cor_1, err_cor_2)
        return noc, insert_header, error_correction

    def remove_channel(self, channel_no):
        """
        To remove a channel from the multiplexer it's exclusive chunks will be allowed for all channels. This prevents
        all chunks to be exclusive by excessive use of add/remove channel.
        :param channel_no:
        :return:
        """
        excl_ch = self.excl_chunks[channel_no]
        self.alld_chunks[channel_no] = set()
        self.excl_chunks[channel_no] = set()
        av_channel = [ind for ind, ch in enumerate(self.channel) if
                      ind != channel_no and len(self.alld_chunks[ind]) > 0]
        if self.secure and self.no_channel < 4:
            print("Removing a channel in secure mode is only possible for 4 or more channel.")
        for ch in self.alld_chunks:
            ch.update(excl_ch)
        print("Removed channel " + str(channel_no) + ".")

    def remove_channel_old(self, channel_no):
        """
        To remove a channel from the multiplexer it's chunk sets will be distributed to the remaining channels.
        :param channel_no:
        :return:
        """
        alld_ch = self.alld_chunks[channel_no]
        excl_ch = self.excl_chunks[channel_no]
        self.alld_chunks[channel_no] = set()
        self.excl_chunks[channel_no] = set()
        av_channel = [ind for ind, ch in enumerate(self.channel) if
                      ind != channel_no and len(self.alld_chunks[ind]) > 0]
        if self.secure and self.no_channel < 4:
            print("Removing a channel in secure mode is only possible for 4 or more channel.")
        spl_excl = self.split_set(excl_ch, self.no_channel - 1)
        for ex in spl_excl:
            ch = random.choice(av_channel)
            self.alld_chunks[ch].update(ex)
            self.excl_chunks[ch].update(ex)
            alld_ch -= ex
        spl_allwd = self.split_set(alld_ch, len(av_channel))
        for allwd in spl_allwd:
            self.alld_chunks[random.choice(av_channel)].update(allwd)
        print("Removed channel " + str(channel_no) + ".")

    def add_channel(self):
        """
        Adds a channel to the multiplexer.
        :return:
        """
        self.no_channel += 1
        self.channel.append(RU10Decoder.pseudo_decoder(self.no_chunks, False))
        self.alld_chunks.append(set())
        self.excl_chunks.append(set())
        if self.secure:
            self.used_chunks.append([])
        chunks = set([x for x in range(0, self.no_chunks)])
        for ex in self.excl_chunks:
            chunks -= ex
        self.excl_chunks[-1].update(set(random.sample(chunks, self.no_excl)))
        for allwd in self.alld_chunks:
            allwd -= self.excl_chunks[-1]
        self.alld_chunks[-1].update(chunks)
        print("Successfully added channel.")

    @staticmethod
    def split_set(dset, no_subsets):
        """
        Splits sets into no_subsets equal sized subsets.
        :param dset:
        :param no_subsets:
        :return:
        """
        res = [set() for _ in range(0, no_subsets)]
        while len(dset) % no_subsets != 0:
            el = random.sample(dset, 1)[0]
            dset.remove(el)
            res[random.randint(0, no_subsets - 1)].add(el)
        slen = int(len(dset) / no_subsets)
        for s in range(0, no_subsets):
            res[s].update(set(random.sample(dset, slen)))
            dset -= res[s]
        return res

    def reset_dec(self):
        """
        Generates a new decoder to reset it.
        :return:
        """
        self.decoder = RU10Decoder.pseudo_decoder(self.encoder.number_of_chunks, False)
        if self.decoder.distribution is None:
            self.decoder.distribution = RaptorDistribution(self.encoder.number_of_chunks)
            self.decoder.number_of_chunks = self.encoder.number_of_chunks
            _, self.decoder.s, self.decoder.h = intermediate_symbols(self.encoder.number_of_chunks,
                                                                     self.decoder.distribution)
            self.decoder.createAuxBlocks()
        print("Resetted the overall decoder.")

    def create_all_packets(self, spare1core=True):
        """
        Creates all possible packets (based on the id_len_format) for the file. Uses multiprocessing.
        :param spare1core:
        :return:
        """
        cores = multiprocessing.cpu_count()
        if spare1core:
            cores -= 1
        p = multiprocessing.Pool(cores)
        max_num = self.max_num_from_id_format(self.id_len_format)
        stepsize = max_num / cores
        param = [int(floor(i * stepsize)) for i in range(cores)]
        while_count = int(ceil(stepsize)) + 1
        a = p.map(partial(self.run, file=self.file, while_count=while_count, number_of_chunks=self.no_chunks,
                          id_len_format=self.id_len_format), param)
        packets = []
        for lst in a:
            packets.extend(lst)
        print("Created " + str(len(packets)) + " packets.")
        return packets

    @staticmethod
    def max_num_from_id_format(id_len_format):
        """
        Calculates the maximal number that can be stored with the id length format.
        :param id_len_format:
        :return:
        """
        return 2 ** (struct.calcsize(id_len_format) * 8)

    @staticmethod
    def run(seq_seed=None, file='Dorn', while_count=1000, number_of_chunks=300, id_len_format='H'):
        """
        Generates packets.
        :param seq_seed:
        :param file:
        :param while_count:
        :param number_of_chunks:
        :param id_len_format:
        :return:
        """
        dist = RaptorDistribution(number_of_chunks)
        x = RU10Encoder(file, number_of_chunks, dist, insert_header=False, id_len_format=id_len_format,
                        number_of_chunks_len_format="B",
                        save_number_of_chunks_in_packet=False)
        x.prepare()
        max_num = Multiplexer.max_num_from_id_format(id_len_format)
        i = 0
        tmp_list = []
        while i < while_count:
            if seq_seed + i >= max_num:
                break
            packet = x.create_new_packet(seed=seq_seed + i)
            tmp_list.append(packet)
            i += 1
            if i % 1000 == 0:
                print(str(i))
        return [ParallelPacket.from_packet(p) for p in tmp_list]


if __name__ == "__main__":
    mltp = Multiplexer('../.INFILES/Dorn', 5, secure=True, no_excl=1, id_len_format='H')
    a = mltp.create_all_packets(True)
    mltp.do_multiplex_with_packets(a, pseudo_decode=False)
