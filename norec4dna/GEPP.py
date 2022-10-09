# Partially Based on jgcastro89 's Code https://gist.github.com/jgcastro89/49090cc69a499a129413597433b9baab
import numpy as np
import typing

from norec4dna.helper import xor_numpy

debug = False


def GEPP(a: np.array, b: np.array):
    return GEPP_intern(a, b)  # np.frombuffer(b, dtype="uint8"))


class GEPP_intern:
    """
    Gaussian elimination with partial pivoting.
    input: A is an n x n np matrix
           b is an n x 1 np array
    output: x is the solution of Ax=b
            with the entries permuted in
            accordance with the pivoting
            done by the algorithm
    post-condition: A and b have been modified.
    :return
    """

    def __init__(self, A: np.array, b: np.array):
        self.A: np.array = A  # input: A is an n x n np matrix
        self.b: np.array = np.fromstring(b, dtype='uint8')  # b is an n x 1 np array
        self.chunk_to_used_packets: np.array = np.identity(max(self.A.shape[0], self.A.shape[1]), dtype=bool)  # inverse part
        while len(self.chunk_to_used_packets) < len(self.A):
            self.chunk_to_used_packets = np.vstack((self.chunk_to_used_packets, np.full((1, len(self.A[1])), False)))
        self.packet_mapping: np.array = np.array([x + 1 for x in range(len(self.A))], dtype=np.uint)
        self.n: int = 0  # n is the length of A
        self.m: int = 0  # m is the width of A
        self.tmp_A: typing.List = []
        self.tmp_b: typing.List = []
        self._update_input()  # method that validates input
        self.result_mapping: np.array = np.zeros((np.int64(self.m), np.int64(1)), np.int32)

    def solve(self, partial=False) -> bool:
        #TODO we might want to remove the added lines if partial=True to avoid unexpected behaviour
        # if the same instance is used with partial=False afterwards
        # additionally the additional lines will most likely slow down the algorithm
        if not (self.isPotentionallySolvable() or partial):
            return False
        try:
            self.insert_tmp()
            # since we need rows >= columns to perform partial pivoting, we fill with full-false rows..
            diff = len(self.A[0]) - len(self.A)
            if diff > 0 and partial:
                for _ in range(diff):
                    self.A = np.vstack((self.A, np.full((1, len(self.A[0])), False, bool)))
                    self.b = np.vstack((self.b, np.array([0 for x in range(len(self.b[0]))], dtype="uint8")))
                    self.packet_mapping = np.concatenate((self.packet_mapping, [len(self.A) + 1]))
                self._update_input()
            return self._elimination()
        except Exception as ex:
            print(f"GEPP: {ex}")
            return False

    def addRow(self, row: typing.Union[typing.Set[int], np.array], data: np.array, ):
        if self.n == 0:
            # first insert
            self.A = np.vstack((self.A, row))
            self.b = np.vstack((self.b, data))
            self.packet_mapping = np.vstack((self.packet_mapping, [len(self.A)]))
            self._update_input()
        self.tmp_A.append(row)
        self.tmp_b.append(data)
        self._update_input()

    def insert_tmp(self):
        if len(self.tmp_A) > 0:
            self.packet_mapping = np.concatenate(
                (self.packet_mapping, np.array([x for x in range(len(self.A) + 1, self.m + len(self.tmp_A) + 1)])))
            self.A = np.vstack((self.A, self.tmp_A))
            self.tmp_A.clear()
        if len(self.tmp_b) > 0:
            self.b = np.vstack((self.b, self.tmp_b))
            self.tmp_b.clear()
        if max(self.A.shape[0], self.A.shape[1]) > self.chunk_to_used_packets.shape[0]:
            # we have to make sure that chunk_to_used_packets.shape is always >= max(A.shape)
            # since we might run solve() multiple times with rows added we can't simply recreate chunks_to_used_packets
            # as a new identity matrix on the fly...
            tmp = np.identity(len(self.A), dtype=bool)
            # create a view overlapping the new matrix with the old one for the top left corner
            tmp[0:self.chunk_to_used_packets.shape[0], 0:self.chunk_to_used_packets.shape[1]] = self.chunk_to_used_packets

            self.chunk_to_used_packets = tmp

    def _update_input(self):
        self.n = len(self.A) + len(self.tmp_A)  # Number of Line
        self.m = len(self.A[0])  # Number of Columns

    def isPotentionallySolvable(self):
        return self.m <= self.n

    def _py_elimination(self) -> bool:
        """
        k represents the current pivot row. Since GE traverses the matrix in the
        upper right triangle, we also use k for indicating the k-th diagonal
        column index.
        :return
        """
        # Create Order Vector
        # self.order = np.array([[i] for i in range(0, self.n)])
        # Elimination
        if self.chunk_to_used_packets is None:
            # TODO OR if len(self.chunk_to_used_packets) < len(self.A) fill with Identity rows
            self.chunk_to_used_packets: np.array = np.identity(max(self.A.shape()[0], self.A.shape()[1]), dtype=np.bool)  # inverse part
            while len(self.chunk_to_used_packets) < self.A.shape()[1]:
                self.chunk_to_used_packets = np.vstack(
                    (self.chunk_to_used_packets, np.full((1, len(self.A[1])), False)))
        for k in range(self.m):  # - 1):
            # Pivot
            maxindex = abs(self.A[k:, k]).argmax() + k
            if self.A[maxindex, k] == 0:
                continue
                # Matrix is singular - raise ValueError("Matrix is singular.")
            # Swap
            if maxindex != k:
                self.A[[k, maxindex]] = self.A[[maxindex, k]]
                self.b[[k, maxindex]] = self.b[[maxindex, k]]
                # swap inverse (this is required if matrix is not square)
                self.chunk_to_used_packets[[k, maxindex]] = self.chunk_to_used_packets[[maxindex, k]]
                tmp = self.packet_mapping[maxindex]
                self.packet_mapping[maxindex] = self.packet_mapping[k]
                self.packet_mapping[k] = tmp
            # Eliminate top to bottom:
            for row in range(k + 1, self.n):
                if self.A[row, k:][0] == 1:
                    self.A[row, k:] = xor_numpy(self.A[row, k:], self.A[k, k:])
                    # Same for data:
                    self.b[row] = xor_numpy(self.b[row], self.b[k])
                    # same for inverse part
                    # TODO: we not not want to go beyond the width of the matrix
                    # (inverse is only defined for square matrices)
                    # ideally we want to only reorder A and b and then eliminate from top to width and from width to top
                    self.chunk_to_used_packets[row] = xor_numpy(self.chunk_to_used_packets[row],
                                                                self.chunk_to_used_packets[k])

        # Eliminate bottom to top:
        for k in range(self.m - 1, -1, -1):
            for row in range(k, -1, -1):
                if self.A[row, k] == 1 and row != k:
                    self.A[row, k:] = xor_numpy(self.A[row, k:], self.A[k, k:])
                    # Same for data:
                    self.b[row] = xor_numpy(self.b[row], self.b[k])
                    # same for inverse part
                    self.chunk_to_used_packets[row] = xor_numpy(self.chunk_to_used_packets[row],
                                                                self.chunk_to_used_packets[k])

        self.result_mapping = self.generateResultMapping()
        return self.isSolved()

    try:
        from cdnarules import elimination  # just to trigger exception before running...

        def _elimination(self) -> bool:
            from cdnarules import elimination
            if self.chunk_to_used_packets is None:
                # TODO OR if len(self.chunk_to_used_packets) < len(self.A) fill with Identity rows
                self.chunk_to_used_packets: np.array = np.identity(max(self.A.shape()[0], self.A.shape()[1]), dtype=bool)  # inverse part
            elimination(self.A, self.b, self.packet_mapping, self.chunk_to_used_packets)
            self.result_mapping = self.generateResultMapping()
            return self.isSolved()
    except:
        print("Gaussian Elimination - C Module failed to load, falling back to slow mode")

        def _elimination(self) -> bool:
            return self._py_elimination()

    def generateResultMapping(self) -> np.array:
        """
        returns which row maps to which raw-data-Chunk. 
        """
        res_rows = np.full((self.m, 1), -1, int)
        for k in range(self.n):
            row = self.A[k]
            if np.sum(row) == 1:
                x = np.where(row == True)[0][0]
                res_rows[x] = k
        return res_rows

    def isSolved(self) -> bool:
        all_true = set(x for x in range(self.m))
        res = len(all_true) == len(set(x[0] for x in self.result_mapping))
        if debug and not res:
            print("Chunk(s) " + str(all_true - set(x[0] for x in self.result_mapping)) + " are missing. ")
        return res

    def getSolvedCount(self) -> int:
        return len(set(x[0] for x in self.result_mapping))

    def get_common_packets(self, chunk_id_lst: typing.List[int], valid_chunks_lst: typing.List[int] = None) -> np.array:
        """
        returns a list of packets that are used to reconstruct all chunks in chunk_id_lst
        :param chunk_id_lst: list of chunk ids that contain invalid data (corrupted)
        :param valid_chunks_lst: list of chunk ids that contain valid data (not corrupted)
        """
        # TODO: if there is no combination of packets that were used for all defined chunks, then return the minimum list of packets that were used for these chunks
        # additionally one could define an additional paramter with a list of all known _GOOD_ chunks,
        # this could then be used as a negative list of packets -> we know that all packets involved in the decoding of that packet must be correct
        if valid_chunks_lst is None:
            valid_chunks_lst = []
        if self.chunk_to_used_packets is None:
            return None
        res = np.ones(len(self.chunk_to_used_packets[0]), dtype=bool)
        for i in chunk_id_lst:
            res = np.logical_and(res, self.chunk_to_used_packets[i])
        for i in valid_chunks_lst:
            res = np.logical_and(res, np.invert(self.chunk_to_used_packets[i]))
        return res

    def remove_row(self, row_id):
        """
        removes the row with the given id from the matrix
        - we might use this option to remove packets that were detected as invalid [manual or (semi-)automatic]
        """
        self.A = np.delete(self.A, row_id, 0)
        self.b = np.delete(self.b, row_id, 0)
        self.packet_mapping = np.delete(self.packet_mapping, row_id, 0)
        self._update_input()
        self._elimination()

    def find_missing_chunks(self):
        """
        returns all chunks from A that do not appear in any received packet and would thus reduce to the complete result
        """
        res = np.zeros(len(self.A[0]), dtype=bool)
        for line in self.A:
            res = np.logical_or(res, line)
        return np.invert(res)


def main():
    A = np.array([[0, 1, 0, 1], [0, 0, 1, 0], [0, 1, 1, 0], [0, 1, 0, 0], [1, 1, 1, 0]], dtype=bool)

    b = np.array(
        [
            [np.frombuffer(b"Hallo", dtype="uint8")],
            [np.frombuffer(b"Welt!", dtype="uint8")],
            [np.frombuffer(b"Test3", dtype="uint8")],
            [np.frombuffer(b"Test1", dtype="uint8")],
            [np.frombuffer(b"12345", dtype="uint8")],
        ]
    )

    gauss_elim_piv = GEPP(np.copy(A), np.copy(b.reshape(len(b), -1)))

    mapping = gauss_elim_piv.generateResultMapping()
    print(mapping)
    gauss_elim_piv.solve()
    # gauss_elim_piv.solve()
    print("Is solved: " + str(gauss_elim_piv.isSolved()))
    print(gauss_elim_piv.A)
    print([i[0] for i in gauss_elim_piv.result_mapping])
    print("Packets used for each chunk:")
    print(gauss_elim_piv.chunk_to_used_packets)
    print("Common packets:")
    print(gauss_elim_piv.get_common_packets([1, 3], [2]))
    print("Missing chunks:")
    print(gauss_elim_piv.find_missing_chunks())
    # print([chr(x) for x in gauss_elim_piv.b[2][0]])

    b_copy = np.copy(gauss_elim_piv.b).reshape(gauss_elim_piv.chunk_to_used_packets.shape[0], -1)
    gauss_elem_inverse = GEPP(np.copy(gauss_elim_piv.chunk_to_used_packets), b_copy)
    gauss_elem_inverse.solve()
    print("Is solved: " + str(gauss_elem_inverse.isSolved()))
    print("".join([chr(x) for x in gauss_elem_inverse.b.reshape(-1)]))

#TODO create a function that (interactivley) shows which chunks could be affected by a given list of packets
# that way we can guide the user to which chunks might also be invalid so that the user can further finetune the search
# and won't waste time tagging packets as valid that are not even in consideration

# TODO: create function that takes a list of chunks and returns the packets that were involoved in the decoding of ALL(!) of these chunks
# this is needed to find out which packet was responsible for the invalid chunks


# TODO: create a function that analyses GEPP.A and return the missing chunks/packets that would allow to solve the system
if __name__ == "__main__":
    main()
