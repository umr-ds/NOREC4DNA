import socket
import struct
import typing


class MultiInterfaceBase:
    __IP_COMMAND = 0x8915
    __BROADCAST_COMMAND = 0x8919
    __NETMASK_COMMAND = 0x891B

    """
    @staticmethod
    def map_ip_interface(ips, interfaces):
        iface_mapping =
        for ip in ips:
            netmask = MultiInterfaceBase.get_ip_address(ip, get_netmask=True)
        return None
    """

    @staticmethod
    def list_interfaces() -> typing.List[str]:
        try:
            return [x[1] for x in socket.if_nameindex()]
        except Exception:
            import netifaces as ni

            return ni.interfaces()

    @staticmethod
    def filter_interfaces(
        ifaces: typing.List[str], get_broadcast: bool = False
    ) -> typing.List[str]:
        clean_ifaces: typing.List[str] = []
        for x in ifaces:
            try:
                print("%s: %s" % (x, MultiInterfaceBase.get_ip_address(x, get_broadcast)))
                if x != "lo":
                    clean_ifaces.append(x)
            except Exception:
                continue
        return clean_ifaces

    @staticmethod
    def get_ip_address(
        ifname: typing.Union[str, bytes], get_broadcast: bool = False, get_netmask: bool = False
    ) -> str:
        assert (get_broadcast is get_netmask and get_netmask is False) or (
            get_netmask is not get_broadcast
        )
        # only broadcast OR netmask OR neither of both are allowed
        ifname_str = ifname.decode() if isinstance(ifname, bytes) else ifname
        try:
            import fcntl

            s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
            cmd = (
                MultiInterfaceBase.__BROADCAST_COMMAND
                if get_broadcast
                else MultiInterfaceBase.__IP_COMMAND
            )
            cmd = MultiInterfaceBase.__NETMASK_COMMAND if get_netmask else cmd
            return socket.inet_ntoa(
                fcntl.ioctl(s.fileno(), cmd, struct.pack("256s", ifname_str[:15].encode()))[20:24]
            )
        except ImportError:
            # non UNIX systems:
            # from netifaces import AF_INET, AF_INET6, AF_LINK, AF_BRIDGE
            import netifaces as ni

            mode = "broadcast" if get_broadcast else "addr"
            mode = "netmask" if get_netmask else mode
            return ni.ifaddresses(ifname_str)[ni.AF_INET][0][mode]
