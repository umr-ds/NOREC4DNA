#!/usr/bin/python
# -*- coding: latin-1 -*-
import binascii
import os
import struct
import sys
import time

import typing
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler

from norec4dna import Decoder
from norec4dna.LTDecoder import LTDecoder
from norec4dna.OnlineDecoder import OnlineDecoder
from norec4dna.RU10Decoder import RU10Decoder


@DeprecationWarning
class Watcher:
    def __init__(self, directory):
        self.DIRECTORY_TO_WATCH = directory
        self.observer = Observer()
        self.use_crc = False
        self.decoders = dict()
        self.decoder_to_filename = dict()

    def get_or_init_decoder(self, name="default") -> typing.Union[Decoder, LTDecoder, RU10Decoder, OnlineDecoder]:
        if name not in self.decoders.keys():
            self.decoder_to_filename[name] = []
            if isinstance(name, int):
                num = int(name)
                # most significant bit indicates which decoder to use:
                tmp = 0b1100000000000000 & num
                if tmp == 0b1000000000000000 or tmp == 0b1100000000000000:
                    # since RU10 looks like the best encoding, we give it most space
                    self.decoders[name] = RU10Decoder()
                elif tmp == 0b0100000000000000:
                    self.decoders[name] = LTDecoder()
                else:  # tmp == 0b0000000000000000:
                    self.decoders[name] = OnlineDecoder()
            self.decoders[name].set_read_all_before_decode(False)
        return self.decoders[name]

    def handle_data(self, data, filename, remove_finished=False):
        dsplit = data.split(",")
        command = dsplit[0]
        if not command.startswith("+RX"):
            return 0
        payload = dsplit[1]

        csplit = command.split(" ")
        command, payload_size = csplit[0], csplit[1]
        if int(payload_size) * 2 != len(payload):
            print(
                "Payload_size={} but actual size is: {}".format(
                    int(payload_size) * 2, len(payload)
                )
            )
            return 0
        rssi, snr = dsplit[2], dsplit[3]
        binary_payload = binascii.unhexlify(payload)
        crc16_name = struct.unpack("!H", binary_payload[:2])[0]
        # print(crc16_name)
        decoder = self.get_or_init_decoder(crc16_name)
        if decoder.GEPP is not None and decoder.is_decoded():
            os.remove(filename)
            try:
                target = decoder.headerChunk.get_file_name().decode()
            except:
                target = ""
            print("{} skipped, target file {} has already been decoded.".format(filename, target))
            return
        print(
            "{} : {} bytes for crc16(filename) = {} - RSSI= {}, SNR={}".format(command, payload_size, crc16_name, rssi,
                                                                               snr))

        self.decoder_to_filename[crc16_name].append(filename)
        packet = decoder.parse_raw_packet(binary_payload[2:], use_crc=self.use_crc)
        if packet is not None:
            decoder.input_new_packet(packet)
        if decoder.number_of_chunks != 0 and decoder.GEPP.n >= decoder.number_of_chunks:
            print(
                "received {} Packets for File consisting of {} chunks.".format(decoder.GEPP.n,
                                                                               decoder.number_of_chunks))
        if decoder.GEPP.isPotentionallySolvable() and decoder.solve():
            decoder.saveDecodedFile()
            if remove_finished:
                for filename in self.decoder_to_filename[crc16_name]:
                    try:
                        os.remove(filename)
                    except:
                        continue
                self.decoder_to_filename[crc16_name].clear()

    @staticmethod
    def load_file(path):
        with open(path, "r") as f:
            content = f.read()
        return content

    def run(self):
        event_handler = Handler()
        event_handler.register("modified", lambda x: self.handle_data(self.load_file(x), x, remove_finished=True), )
        self.observer.schedule(event_handler, self.DIRECTORY_TO_WATCH, recursive=False)
        self.observer.start()
        try:
            while True:
                time.sleep(5)
        except:
            self.observer.stop()
            print("Error")
        self.observer.join()


class Handler(FileSystemEventHandler):
    def __init__(self):
        self.registry = {"created": [], "modified": [], "moved": [], "deleted": []}

    def on_any_event(self, event):
        if event.is_directory:
            return None
        for function in self.registry[event.event_type]:
            function(event.src_path)

    def register(self, event_type, function):
        """
        event_type might be "modified", "created", "moved" or "deleted" (or "all")
        function should have exactly one parameter: f(event.src_path)
        """
        if event_type == "all":
            for key in self.registry.keys():
                self.registry[key].append(function)
        else:
            self.registry[event_type].append(function)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python {} <folder to watch>".format(sys.argv[0].split("/")[-1]))
        exit()
    w = Watcher(sys.argv[1])
    w.run()
