import socket
import struct
import threading
from queue import Empty, Queue

from ..ErrorCorrection import reed_solomon_decode
from ..helper import xor_mask
from ..RU10Decoder import RU10Decoder
from .MultiInterfaceBase import MultiInterfaceBase

ANY_INTERFACE_IP = socket.inet_ntoa(struct.pack("!I", socket.INADDR_ANY))


class MultiInterfaceReceiver(MultiInterfaceBase):
    def __init__(self, listen_ifaces=None, port=45155):
        self.listen_ifaces = listen_ifaces
        self.__port = port

    def get_port(self):
        return self.__port

    def create_listen_socket(self, interface, broadcast=False):
        sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        ip = self.get_ip_address(interface, broadcast)
        if ip == ANY_INTERFACE_IP:
            return None
        sock.bind((ip, self.__port))
        print("Listening on %s:%s" % (ip, self.get_port()))
        # sock.setblocking(False)
        sock.settimeout(1)
        return sock


def listen(sock, queue, signals):
    if sock is None:
        return None
    while not signals["shutdown"]:
        try:
            data, sender = sock.recvfrom(1024)
            a, b = struct.unpack("<II", data[0:8])
            print(
                "Packet from: %s:%s to %s - #Chunks: %s - Id: %s"
                % (sender[0], sender[1], sock.getsockname()[0], xor_mask(a), xor_mask(b))
            )
            queue.put(data)
        except socket.error:
            # queue.put(e)
            continue


def _create_sockets(receiver, interfaces, broadcast):
    socks = []
    for interface in interfaces:
        try:
            socks.append(receiver.create_listen_socket(interface, broadcast))
        except Exception as exc:
            print("<%s>: %s" % (interface, exc))
            raise exc
    return socks


def _start_listener_threads(socks, pqueue, signals) -> None:
    for sock in socks:
        try:
            thread = threading.Thread(target=listen, args=(sock, pqueue, signals))
            thread.start()
        except socket.error as exc:
            print(exc)


def _get_packet_batch(pqueue, sock_count):
    packet_strs = []
    while True:
        try:
            packet_str = pqueue.get(timeout=2)
            packet_strs.append(packet_str)
            if len(packet_strs) > 50 * sock_count:
                return packet_strs
        except Empty as exc:
            print(exc)
            return packet_strs


def _process_packet_batch(decoder, packet_strs):
    for packet_str in packet_strs:
        pack = decoder.parse_raw_packet(
            packet_str,
            crc_len_format="L",
            number_of_chunks_len_format="I",
            packet_len_format="I",
            id_len_format="I",
        )
        if isinstance(pack, str):
            continue
        decoder.input_new_packet(pack)


def _run_receiver() -> None:
    receiver = MultiInterfaceReceiver()
    broadcast = False
    use_header_chunk = True
    decoder = RU10Decoder(
        None,
        use_headerchunk=use_header_chunk,
        error_correction=reed_solomon_decode,
        static_number_of_chunks=None,
    )
    decoder.read_all_before_decode = True
    clean_ifaces = receiver.filter_interfaces(receiver.list_interfaces())
    socks = _create_sockets(receiver, clean_ifaces, broadcast)
    pqueue = Queue()
    signals = {"shutdown": False}
    _start_listener_threads(socks, pqueue, signals)
    while True:
        try:
            packet_strs = _get_packet_batch(pqueue, len(socks))
            _process_packet_batch(decoder, packet_strs)
            if decoder.GEPP is not None and decoder.solve():
                print("Success!")
                signals["shutdown"] = True
                pqueue.task_done()
                decoder.saveDecodedFile(null_is_terminator=False, print_to_output=False)
                break
        except Exception as exc:
            print(exc)


if __name__ == "__main__":
    _run_receiver()
