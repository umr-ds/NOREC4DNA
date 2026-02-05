#if defined(__AVX2__)
#include <immintrin.h>
#define USE_AVX2 1
#elif defined(__SSE2__)
#include <emmintrin.h>
#define USE_SSE2 1
#endif

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <limits.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

// Cross-compiler portability shims (MSVC vs GCC/Clang)
#if defined(_MSC_VER)
  #include <intrin.h>    // __popcnt64
  #include <malloc.h>    // _aligned_malloc/_aligned_free
  #ifndef ALIGN_BYTES
  #define ALIGN_BYTES 32
  #endif
  #define ALIGNED __declspec(align(ALIGN_BYTES))
  #define RESTRICT __restrict
  #define INLINE static __forceinline
  #define LIKELY(x)   (x)
  #define UNLIKELY(x) (x)
  static void* aligned_malloc(size_t size) { return _aligned_malloc(size, ALIGN_BYTES); }
  static void aligned_free(void* p) { _aligned_free(p); }
#else
  #ifndef ALIGN_BYTES
  #define ALIGN_BYTES 32
  #endif
  #define ALIGNED __attribute__((aligned(ALIGN_BYTES)))
  #define RESTRICT __restrict__
  #define INLINE static inline __attribute__((always_inline))
  #define LIKELY(x)   __builtin_expect(!!(x), 1)
  #define UNLIKELY(x) __builtin_expect(!!(x), 0)
  static void* aligned_malloc(size_t size) {
      void* ptr = NULL;
      if (posix_memalign(&ptr, ALIGN_BYTES, size) != 0) return NULL;
      return ptr;
  }
  static void aligned_free(void* p) { free(p); }
#endif

// Fast array pointer access with stride calculation
INLINE void* fast_array_ptr2(PyArrayObject* arr, npy_intp i, npy_intp j) {
    char* data = (char*)PyArray_DATA(arr);
    npy_intp stride0 = PyArray_STRIDE(arr, 0);
    npy_intp stride1 = PyArray_STRIDE(arr, 1);
    return (void*)(data + i * stride0 + j * stride1);
}

// Type aliases used throughout
#define T uint64_t
#define BYTE uint8_t


INLINE void* fast_array_ptr1(PyArrayObject* arr, npy_intp i) {
    char* data = (char*)PyArray_DATA(arr);
    npy_intp stride0 = PyArray_STRIDE(arr, 0);
    return (void*)(data + i * stride0);
}

static PyObject* bitSet(PyObject* self, PyObject *args) {
    T v;
    unsigned int b;
    if (UNLIKELY(!PyArg_ParseTuple(args, "Ki", &v, &b))) {
        return NULL;
    }
    if (UNLIKELY(b >= 64)) {
        PyErr_SetString(PyExc_ValueError, "Bit position must be less than 64");
        return NULL;
    }
    bool c = ((v >> b) & 1) == 1;
    return PyBool_FromLong(c);
}

// Optimized bit counting using builtin when available
INLINE int bitsSet_internal(T v) {
#if defined(_MSC_VER)
    return (int)__popcnt64((unsigned __int64)v);
#elif defined(__BUILTIN_POPCOUNT)
    return __builtin_popcountll(v);
#else
    T x = v;
    x = x - ((x >> 1) & (T)~(T)0/3);
    x = (x & (T)~(T)0/15*3) + ((x >> 2) & (T)~(T)0/15*3);
    x = (x + (x >> 4)) & (T)~(T)0/255*15;
    return (T)(x * ((T)~(T)0/255)) >> (sizeof(T) - 1) * CHAR_BIT;
#endif
}


static PyObject* bitsSet(PyObject* self, PyObject *args) {
    T v;
    if (UNLIKELY(!PyArg_ParseTuple(args, "K", &v))) {
        return NULL;
    }
    int c = bitsSet_internal(v);
    return PyLong_FromLong(c);
}

INLINE T grayCode_internal(T x) {
    return (x >> 1) ^ x;
}

static PyObject* grayCode(PyObject* self, PyObject *args) {
    uint64_t x;
    if (UNLIKELY(!PyArg_ParseTuple(args, "K", &x))) {
        return NULL;
    }
    T c = grayCode_internal(x);
    return PyLong_FromUnsignedLongLong(c);
}

static PyObject* buildGraySequence(PyObject* self, PyObject *args) {
    int length, b;
    if (UNLIKELY(!PyArg_ParseTuple(args, "ii", &length, &b))) {
        return NULL;
    }
    if (UNLIKELY(length <= 0)) {
        PyErr_SetString(PyExc_ValueError, "Length must be positive");
        return NULL;
    }

    npy_intp dims = length;
    PyObject *result = PyArray_SimpleNew(1, &dims, NPY_ULONGLONG);
    if (UNLIKELY(!result)) {
        return PyErr_NoMemory();
    }

    T* RESTRICT resultDataPtr = (T*)PyArray_DATA((PyArrayObject*)result);
    T x = 0;
    int i = 0;

    while (i < length) {
        T g = grayCode_internal(x);
        if (LIKELY(bitsSet_internal(g) == b)) {
            resultDataPtr[i++] = g;
        }
        x++;
        // Prevent infinite loop for impossible combinations
        if (UNLIKELY(x == 0)) { // Overflow check
            break;
        }
    }
    return result;
}

// Optimized XOR operations with vectorization hints
static void do_xor_bool_optimized(bool* RESTRICT arr_a, bool* RESTRICT arr_b,
                                 npy_intp length, bool* RESTRICT outArr) {
    npy_intp i;
    #if defined(__GNUC__)
    #pragma GCC ivdep
    #pragma GCC vector
    #endif // __GNUC__
    for (i = 0; i < length; i++) {
        outArr[i] = arr_a[i] ^ arr_b[i];
    }
}

static void do_xor_byte_optimized(BYTE* RESTRICT arr_a, BYTE* RESTRICT arr_b,
                                 npy_intp length, BYTE* RESTRICT outArr) {
#if defined(USE_AVX2)
    npy_intp i = 0;
    // Process 32 bytes at a time with AVX2
    npy_intp avx_end = (length / 32) * 32;
    for (i = 0; i < avx_end; i += 32) {
        __m256i a = _mm256_loadu_si256((__m256i*)(arr_a + i));
        __m256i b = _mm256_loadu_si256((__m256i*)(arr_b + i));
        __m256i result = _mm256_xor_si256(a, b);
        _mm256_storeu_si256((__m256i*)(outArr + i), result);
    }
    // Handle remaining bytes
    for (; i < length; i++) {
        outArr[i] = arr_a[i] ^ arr_b[i];
    }
#elif defined(USE_SSE2)
    npy_intp i = 0;
    // Process 16 bytes at a time with SSE2
    npy_intp sse_end = (length / 16) * 16;
    for (i = 0; i < sse_end; i += 16) {
        __m128i a = _mm_loadu_si128((__m128i*)(arr_a + i));
        __m128i b = _mm_loadu_si128((__m128i*)(arr_b + i));
        __m128i result = _mm_xor_si128(a, b);
        _mm_storeu_si128((__m128i*)(outArr + i), result);
    }
    // Handle remaining bytes
    for (; i < length; i++) {
        outArr[i] = arr_a[i] ^ arr_b[i];
    }
#else
    npy_intp i;
    // Fallback to 64-bit chunks
    npy_intp chunk_size = length & ~7;
    for (i = 0; i < chunk_size; i += 8) {
        uint64_t* a64 = (uint64_t*)(arr_a + i);
        uint64_t* b64 = (uint64_t*)(arr_b + i);
        uint64_t* out64 = (uint64_t*)(outArr + i);
        *out64 = *a64 ^ *b64;
    }
    for (i = chunk_size; i < length; i++) {
        outArr[i] = arr_a[i] ^ arr_b[i];
    }
#endif
}

static void printBin(BYTE in) {
    for (int i = 0; i < 8; i++) {
        PySys_WriteStdout("%d ", (bool)(in & (1 << i)));
    }
}

static void byte_printmatrix(PyArrayObject* A) {
    npy_intp dims_a_0 = PyArray_DIM(A, 0);
    npy_intp dims_a_1 = PyArray_DIM(A, 1);
    for (npy_intp i = 0; i < dims_a_0; i++) {
        PySys_WriteStdout("| ");
        for (npy_intp j = 0; j < dims_a_1; j++) {
            PySys_WriteStdout("%i ", *((BYTE*)fast_array_ptr2(A, i, j)));
        }
        PySys_WriteStdout(" |\n");
    }
}

static void printmatrix(PyArrayObject* A) {
    npy_intp dims_a_0 = PyArray_DIM(A, 0);
    npy_intp dims_a_1 = PyArray_DIM(A, 1);
    for (npy_intp i = 0; i < dims_a_0; i++) {
        PySys_WriteStdout("| ");
        for (npy_intp j = 0; j < dims_a_1; j++) {
            PySys_WriteStdout("%d ", *((bool*)fast_array_ptr2(A, i, j)));
        }
        PySys_WriteStdout(" |\n");
    }
}

INLINE bool isSolvable(PyArrayObject* A) {
    return PyArray_DIM(A, 0) >= PyArray_DIM(A, 1);
}

static PyObject* elimination(PyObject *self, PyObject *args) {
    PyArrayObject *A, *b, *packet_mapping, *chunk_to_used_packets;
    if (UNLIKELY(!PyArg_ParseTuple(args, "O!O!O!O!", &PyArray_Type, &A, &PyArray_Type, &b,
                                  &PyArray_Type, &packet_mapping, &PyArray_Type, &chunk_to_used_packets))) {
        return NULL;
    }

    if (UNLIKELY(!isSolvable(A))) {
        Py_RETURN_FALSE;
    }

    npy_intp dims_a_0 = PyArray_DIM(A, 0); // rows
    npy_intp dims_a_1 = PyArray_DIM(A, 1); // columns
    npy_intp dims_b_1 = PyArray_DIM(b, 1);
    npy_intp dims_chunk_to_used_packets_1 = PyArray_DIM(chunk_to_used_packets, 1);

    bool *dirty_rows = (bool*)aligned_malloc(dims_a_0 * sizeof(bool));
    if (UNLIKELY(!dirty_rows)) {
        return PyErr_NoMemory();
    }

    // Initialize dirty_rows with vectorization hint
    #if defined(__GNUC__)
    #pragma GCC ivdep
    #endif // __GNUC__
    for (npy_intp i = 0; i < dims_a_0; i++) {
        dirty_rows[i] = false;
    }

    bool dirty = false;
    uint8_t num_dirty_rows = 0;

    // Forward elimination
    for (npy_intp i = 0; i < dims_a_1; i++) {
        // Find pivot
        npy_intp pivot = -1;
        for (npy_intp j = i; j < dims_a_0; j++) {
            if (*((bool*)fast_array_ptr2(A, j, i))) {
                pivot = j;
                break;
            }
        }

        if (pivot == -1) {
            PySys_WriteStdout("Could not decode Chunk %" NPY_INTP_FMT "\n", i);
            dirty_rows[i] = true;
            dirty = true;
            num_dirty_rows++;
            continue;
        }

        // Swap rows if needed
        if (pivot != i) {
            // XOR-based row swapping for A
            do_xor_bool_optimized((bool*)fast_array_ptr2(A, i, 0),
                                 (bool*)fast_array_ptr2(A, pivot, 0),
                                 dims_a_1, (bool*)fast_array_ptr2(A, i, 0));
            do_xor_bool_optimized((bool*)fast_array_ptr2(A, i, 0),
                                 (bool*)fast_array_ptr2(A, pivot, 0),
                                 dims_a_1, (bool*)fast_array_ptr2(A, pivot, 0));
            do_xor_bool_optimized((bool*)fast_array_ptr2(A, i, 0),
                                 (bool*)fast_array_ptr2(A, pivot, 0),
                                 dims_a_1, (bool*)fast_array_ptr2(A, i, 0));

            // XOR-based row swapping for b
            do_xor_byte_optimized((BYTE*)fast_array_ptr2(b, i, 0),
                                 (BYTE*)fast_array_ptr2(b, pivot, 0),
                                 dims_b_1, (BYTE*)fast_array_ptr2(b, i, 0));
            do_xor_byte_optimized((BYTE*)fast_array_ptr2(b, pivot, 0),
                                 (BYTE*)fast_array_ptr2(b, i, 0),
                                 dims_b_1, (BYTE*)fast_array_ptr2(b, pivot, 0));
            do_xor_byte_optimized((BYTE*)fast_array_ptr2(b, i, 0),
                                 (BYTE*)fast_array_ptr2(b, pivot, 0),
                                 dims_b_1, (BYTE*)fast_array_ptr2(b, i, 0));

            // XOR-based row swapping for chunk_to_used_packets
            do_xor_bool_optimized((bool*)fast_array_ptr2(chunk_to_used_packets, i, 0),
                                 (bool*)fast_array_ptr2(chunk_to_used_packets, pivot, 0),
                                 dims_chunk_to_used_packets_1,
                                 (bool*)fast_array_ptr2(chunk_to_used_packets, i, 0));
            do_xor_bool_optimized((bool*)fast_array_ptr2(chunk_to_used_packets, i, 0),
                                 (bool*)fast_array_ptr2(chunk_to_used_packets, pivot, 0),
                                 dims_chunk_to_used_packets_1,
                                 (bool*)fast_array_ptr2(chunk_to_used_packets, pivot, 0));
            do_xor_bool_optimized((bool*)fast_array_ptr2(chunk_to_used_packets, i, 0),
                                 (bool*)fast_array_ptr2(chunk_to_used_packets, pivot, 0),
                                 dims_chunk_to_used_packets_1,
                                 (bool*)fast_array_ptr2(chunk_to_used_packets, i, 0));

            // Swap packet_mapping
            PyObject *old_i = PyArray_GETITEM(packet_mapping, (char*)fast_array_ptr1(packet_mapping, i));
            PyObject *old_j = PyArray_GETITEM(packet_mapping, (char*)fast_array_ptr1(packet_mapping, pivot));
            PyArray_SETITEM(packet_mapping, (char*)fast_array_ptr1(packet_mapping, pivot), old_i);
            PyArray_SETITEM(packet_mapping, (char*)fast_array_ptr1(packet_mapping, i), old_j);
            Py_DECREF(old_i);
            Py_DECREF(old_j);
        }

        // Eliminate below pivot
        for (npy_intp j = i + 1; j < dims_a_0; j++) {
            if (dirty_rows[j]) continue;

            if (*((bool*)fast_array_ptr2(A, j, i))) {
                do_xor_bool_optimized((bool*)fast_array_ptr2(A, j, 0),
                                     (bool*)fast_array_ptr2(A, i, 0),
                                     dims_a_1, (bool*)fast_array_ptr2(A, j, 0));
                do_xor_byte_optimized((BYTE*)fast_array_ptr2(b, j, 0),
                                     (BYTE*)fast_array_ptr2(b, i, 0),
                                     dims_b_1, (BYTE*)fast_array_ptr2(b, j, 0));
                do_xor_bool_optimized((bool*)fast_array_ptr2(chunk_to_used_packets, j, 0),
                                     (bool*)fast_array_ptr2(chunk_to_used_packets, i, 0),
                                     dims_chunk_to_used_packets_1,
                                     (bool*)fast_array_ptr2(chunk_to_used_packets, j, 0));
            }
        }
    }

    // Backward elimination
    for (npy_intp col = dims_a_1 - 1; col >= 0; col--) {
        if (dirty_rows[col]) continue;

        for (npy_intp row = col - 1; row >= 0; row--) {
            if (dirty_rows[row]) continue;

            if (*((bool*)fast_array_ptr2(A, row, col))) {
                do_xor_bool_optimized((bool*)fast_array_ptr2(A, row, 0),
                                     (bool*)fast_array_ptr2(A, col, 0),
                                     dims_a_1, (bool*)fast_array_ptr2(A, row, 0));
                do_xor_byte_optimized((BYTE*)fast_array_ptr2(b, row, 0),
                                     (BYTE*)fast_array_ptr2(b, col, 0),
                                     dims_b_1, (BYTE*)fast_array_ptr2(b, row, 0));
                do_xor_bool_optimized((bool*)fast_array_ptr2(chunk_to_used_packets, row, 0),
                                     (bool*)fast_array_ptr2(chunk_to_used_packets, col, 0),
                                     dims_chunk_to_used_packets_1,
                                     (bool*)fast_array_ptr2(chunk_to_used_packets, row, 0));
            }
        }
    }
    aligned_free(dirty_rows);
    return PyBool_FromLong(!dirty);
}

static PyObject* xor_array(PyObject *self, PyObject *args) {
    PyArrayObject *X, *Y;

    if (UNLIKELY(!PyArg_ParseTuple(args, "O!O!", &PyArray_Type, &X, &PyArray_Type, &Y))) {
        return NULL;
    }

    npy_intp dims_x = PyArray_DIM(X, 0);
    if (UNLIKELY(dims_x != PyArray_DIM(Y, 0))) {
        PyErr_SetString(PyExc_ValueError, "input dimensions differ");
        return NULL;
    }

    BYTE* x_DataPtr = (BYTE*)PyArray_DATA(X);
    BYTE* y_DataPtr = (BYTE*)PyArray_DATA(Y);

    npy_intp dims[1] = {dims_x};
    PyObject *out = PyArray_SimpleNew(1, dims, NPY_UINT8);
    if (UNLIKELY(!out)) {
        return PyErr_NoMemory();
    }

    BYTE* out_DataPtr = (BYTE*)PyArray_DATA((PyArrayObject*)out);
    do_xor_byte_optimized(x_DataPtr, y_DataPtr, dims_x, out_DataPtr);

    return out;
}

static PyObject* microsatellite(PyObject* self, PyObject *args) {
    char *text;
    int lengthToLookFor;
    if (UNLIKELY(!PyArg_ParseTuple(args, "si", &text, &lengthToLookFor))) {
        return NULL;
    }
    if (UNLIKELY(lengthToLookFor <= 0)) {
        PyErr_SetString(PyExc_ValueError, "Length must be positive");
        return NULL;
    }

    size_t n = strlen(text);
    if (UNLIKELY(n < 2 * lengthToLookFor)) {
        return Py_BuildValue("(is)", 1, "");
    }

    int res = 1;
    char *resChars = (char*)PyMem_Malloc((lengthToLookFor + 1) * sizeof(char));
    if (UNLIKELY(!resChars)) {
        return PyErr_NoMemory();
    }

    strncpy(resChars, text, lengthToLookFor);
    resChars[lengthToLookFor] = '\0';

    int maxLength = 0;
    size_t i = 0;

    while (i <= n - 2 * lengthToLookFor) {
        if (memcmp(&text[i], &text[i + lengthToLookFor], lengthToLookFor) == 0) {
            res++;
        } else {
            if (maxLength < res) {
                maxLength = res;
                strncpy(resChars, &text[i], lengthToLookFor);
                resChars[lengthToLookFor] = '\0';
            }
            res = 1;
        }
        i += lengthToLookFor;
    }

    if (maxLength < res) {
        maxLength = res;
        strncpy(resChars, &text[i], lengthToLookFor);
        resChars[lengthToLookFor] = '\0';
    }

    PyObject *return_val = Py_BuildValue("(is)", maxLength, resChars);
    PyMem_Free(resChars);
    return return_val;
}

static PyObject* longestSequenceOfChar(PyObject* self, PyObject *args) {
    char *text;
    char *char_x;
    if (UNLIKELY(!PyArg_ParseTuple(args, "ss", &text, &char_x))) {
        return NULL;
    }
    if (UNLIKELY(strlen(char_x) != 1)) {
        PyErr_SetString(PyExc_ValueError, "Second argument must be a single character");
        return NULL;
    }

    int c = 0;
    char res = char_x[0];
    size_t n = strlen(text);
    int curr = 1;

    for (size_t i = 0; i < n - 1; i++) {
        if (text[i] == text[i + 1]) {
            curr++;
        } else {
            if (curr > c && (char_x[0] == '*' || text[i] == char_x[0])) {
                c = curr;
                res = text[i];
            }
            curr = 1;
        }
    }

    if (n > 0 && curr > c && (char_x[0] == '*' || text[n-1] == char_x[0])) {
        c = curr;
        res = text[n-1];
    }

    return Py_BuildValue("(ci)", res, c);
}

static PyObject* repeatRegion(PyObject* self, PyObject *args) {
    char *text;
    int lengthToLookFor;
    if (UNLIKELY(!PyArg_ParseTuple(args, "si", &text, &lengthToLookFor))) {
        return NULL;
    }
    if (UNLIKELY(lengthToLookFor <= 0)) {
        PyErr_SetString(PyExc_ValueError, "Length must be positive");
        return NULL;
    }

    size_t len = strlen(text);
    if (len < lengthToLookFor) {
        return PyLong_FromLong(0);
    }

    char *subseq = (char*)PyMem_Malloc((lengthToLookFor + 1) * sizeof(char));
    if (UNLIKELY(!subseq)) {
        return PyErr_NoMemory();
    }

    int res = 0;
    for (size_t i = 0; i <= len - lengthToLookFor; i++) {
        strncpy(subseq, &text[i], lengthToLookFor);
        subseq[lengthToLookFor] = '\0';
        if (strstr(&text[i + 1], subseq) != NULL) {
            res = 1;
            break;
        }
    }

    PyMem_Free(subseq);
    return PyLong_FromLong(res);
}

static PyObject* smallRepeatRegion(PyObject* self, PyObject *args) {
    char *text;
    int lengthToLookFor;
    if (UNLIKELY(!PyArg_ParseTuple(args, "si", &text, &lengthToLookFor))) {
        return NULL;
    }
    if (UNLIKELY(lengthToLookFor <= 0)) {
        PyErr_SetString(PyExc_ValueError, "Length must be positive");
        return NULL;
    }

    size_t len = strlen(text);
    if (len == 0) {
        return PyFloat_FromDouble(1.0);
    }

    char *subseq = (char*)PyMem_Malloc((lengthToLookFor + 1) * sizeof(char));
    if (UNLIKELY(!subseq)) {
        return PyErr_NoMemory();
    }

    double res = 1.0;
    for (size_t i = 0; i <= len - lengthToLookFor; i++) {
        strncpy(subseq, &text[i], lengthToLookFor);
        subseq[lengthToLookFor] = '\0';
        if (strstr(&text[i + 1], subseq) != NULL) {
            res += 1.0;
        }
    }

    double ratio = res * lengthToLookFor / len;
    if (ratio > 0.44) {
        res = 1.0;
    } else {
        res = ratio * 0.5;
    }

    PyMem_Free(subseq);
    return PyFloat_FromDouble(res);
}

static PyObject* getQUAT(PyObject* self, PyObject *args) {
    int bit1, bit2;
    if (UNLIKELY(!PyArg_ParseTuple(args, "pp", &bit1, &bit2))) {
        return NULL;
    }

    char res;
    if (bit1 && bit2) {
        res = 'T';
    } else if (!bit1 && bit2) {
        res = 'C';
    } else if (bit1 && !bit2) {
        res = 'G';
    } else {
        res = 'A';
    }

    return PyUnicode_FromStringAndSize(&res, 1);
}

static PyObject* byte2QUATS(PyObject* self, PyObject *args) {
    int byte;
    if (UNLIKELY(!PyArg_ParseTuple(args, "i", &byte))) {
        return NULL;
    }

    char res[5];
    int pos = 6;

    for (int i = 0; i < 4; i++) {
        int bit1 = (byte >> (pos + 1)) & 0x01;
        int bit2 = (byte >> pos) & 0x01;
        switch ((bit1 << 1) | bit2) {
            case 0b11: res[i] = 'T'; break;
            case 0b01: res[i] = 'C'; break;
            case 0b10: res[i] = 'G'; break;
            default:   res[i] = 'A'; break;
        }
        pos -= 2;
    }
    res[4] = '\0';

    return PyUnicode_FromString(res);
}

static PyObject* gc_content(PyObject* self, PyObject* args) {
    const char* text;
    if (UNLIKELY(!PyArg_ParseTuple(args, "s", &text))) {
        return NULL;
    }

    size_t count_gc = 0;
    size_t text_length = strlen(text);

    // Vectorization hint for simple counting loop
    #if defined(__GNUC__)
    #pragma GCC ivdep
    #endif // __GNUC__
    for (size_t i = 0; i < text_length; i++) {
        if (text[i] == 'G' || text[i] == 'C') {
            count_gc++;
        }
    }

    if (text_length == 0) {
        return PyFloat_FromDouble(0.0);
    }

    double percentage = (count_gc * 100.0) / text_length;
    return PyFloat_FromDouble(percentage);
}

static PyObject* strContainsSub(PyObject* self, PyObject *args) {
    char *text;
    char *substr;
    if (UNLIKELY(!PyArg_ParseTuple(args, "ss", &text, &substr))) {
        return NULL;
    }

    bool found = strstr(text, substr) != NULL;
    return PyBool_FromLong(found);
}

static PyObject* elimination_with_first_row(PyObject *self, PyObject *args) {
    PyArrayObject *A, *b, *packet_mapping, *chunk_to_used_packets;
    npy_intp first_row_idx = -1;

    if (UNLIKELY(!PyArg_ParseTuple(args, "O!O!O!O!|n", &PyArray_Type, &A, &PyArray_Type, &b,
                                  &PyArray_Type, &packet_mapping, &PyArray_Type, &chunk_to_used_packets,
                                  &first_row_idx))) {
        return NULL;
    }

    if (UNLIKELY(!isSolvable(A))) {
        Py_RETURN_FALSE;
    }

    npy_intp dims_a_0 = PyArray_DIM(A, 0); // rows
    npy_intp dims_a_1 = PyArray_DIM(A, 1); // columns
    npy_intp dims_b_1 = PyArray_DIM(b, 1);
    npy_intp dims_chunk_to_used_packets_1 = PyArray_DIM(chunk_to_used_packets, 1);

    bool *dirty_rows = (bool*)aligned_malloc(dims_a_0 * sizeof(bool));
    if (UNLIKELY(!dirty_rows)) {
        return PyErr_NoMemory();
    }

    // Initialize dirty_rows with vectorization hint
    #pragma GCC ivdep
    for (npy_intp i = 0; i < dims_a_0; i++) {
        dirty_rows[i] = false;
    }

    bool dirty = false;
    uint8_t num_dirty_rows = 0;

    // Forward elimination
    for (npy_intp i = 0; i < dims_a_1; i++) {
        // Find pivot
        npy_intp pivot = -1;

        // If first_row_idx is specified and we're at column 0, use it if valid
        if (i == 0 && first_row_idx >= 0 && first_row_idx < dims_a_0) {
            if (*((bool*)fast_array_ptr2(A, first_row_idx, 0))) {
                pivot = first_row_idx;
            }
        }

        // If no valid first_row_idx or not at column 0, find pivot normally
        if (pivot == -1) {
            for (npy_intp j = i; j < dims_a_0; j++) {
                if (*((bool*)fast_array_ptr2(A, j, i))) {
                    pivot = j;
                    break;
                }
            }
        }

        if (pivot == -1) {
            PySys_WriteStdout("Could not decode Chunk %ld\n", i);
            dirty_rows[i] = true;
            dirty = true;
            num_dirty_rows++;
            continue;
        }

        // Swap rows if needed
        if (pivot != i) {
            // XOR-based row swapping for A
            do_xor_bool_optimized((bool*)fast_array_ptr2(A, i, 0),
                                 (bool*)fast_array_ptr2(A, pivot, 0),
                                 dims_a_1, (bool*)fast_array_ptr2(A, i, 0));
            do_xor_bool_optimized((bool*)fast_array_ptr2(A, i, 0),
                                 (bool*)fast_array_ptr2(A, pivot, 0),
                                 dims_a_1, (bool*)fast_array_ptr2(A, pivot, 0));
            do_xor_bool_optimized((bool*)fast_array_ptr2(A, i, 0),
                                 (bool*)fast_array_ptr2(A, pivot, 0),
                                 dims_a_1, (bool*)fast_array_ptr2(A, i, 0));

            // XOR-based row swapping for b
            do_xor_byte_optimized((BYTE*)fast_array_ptr2(b, i, 0),
                                 (BYTE*)fast_array_ptr2(b, pivot, 0),
                                 dims_b_1, (BYTE*)fast_array_ptr2(b, i, 0));
            do_xor_byte_optimized((BYTE*)fast_array_ptr2(b, pivot, 0),
                                 (BYTE*)fast_array_ptr2(b, i, 0),
                                 dims_b_1, (BYTE*)fast_array_ptr2(b, pivot, 0));
            do_xor_byte_optimized((BYTE*)fast_array_ptr2(b, i, 0),
                                 (BYTE*)fast_array_ptr2(b, pivot, 0),
                                 dims_b_1, (BYTE*)fast_array_ptr2(b, i, 0));

            // XOR-based row swapping for chunk_to_used_packets
            do_xor_bool_optimized((bool*)fast_array_ptr2(chunk_to_used_packets, i, 0),
                                 (bool*)fast_array_ptr2(chunk_to_used_packets, pivot, 0),
                                 dims_chunk_to_used_packets_1,
                                 (bool*)fast_array_ptr2(chunk_to_used_packets, i, 0));
            do_xor_bool_optimized((bool*)fast_array_ptr2(chunk_to_used_packets, i, 0),
                                 (bool*)fast_array_ptr2(chunk_to_used_packets, pivot, 0),
                                 dims_chunk_to_used_packets_1,
                                 (bool*)fast_array_ptr2(chunk_to_used_packets, pivot, 0));
            do_xor_bool_optimized((bool*)fast_array_ptr2(chunk_to_used_packets, i, 0),
                                 (bool*)fast_array_ptr2(chunk_to_used_packets, pivot, 0),
                                 dims_chunk_to_used_packets_1,
                                 (bool*)fast_array_ptr2(chunk_to_used_packets, i, 0));

            // Swap packet_mapping
            PyObject *old_i = PyArray_GETITEM(packet_mapping, fast_array_ptr1(packet_mapping, i));
            PyObject *old_j = PyArray_GETITEM(packet_mapping, fast_array_ptr1(packet_mapping, pivot));
            PyArray_SETITEM(packet_mapping, fast_array_ptr1(packet_mapping, pivot), old_i);
            PyArray_SETITEM(packet_mapping, fast_array_ptr1(packet_mapping, i), old_j);
            Py_DECREF(old_i);
            Py_DECREF(old_j);
        }

        // Eliminate below pivot
        for (npy_intp j = i + 1; j < dims_a_0; j++) {
            if (dirty_rows[j]) continue;

            if (*((bool*)fast_array_ptr2(A, j, i))) {
                do_xor_bool_optimized((bool*)fast_array_ptr2(A, j, 0),
                                     (bool*)fast_array_ptr2(A, i, 0),
                                     dims_a_1, (bool*)fast_array_ptr2(A, j, 0));
                do_xor_byte_optimized((BYTE*)fast_array_ptr2(b, j, 0),
                                     (BYTE*)fast_array_ptr2(b, i, 0),
                                     dims_b_1, (BYTE*)fast_array_ptr2(b, j, 0));
                do_xor_bool_optimized((bool*)fast_array_ptr2(chunk_to_used_packets, j, 0),
                                     (bool*)fast_array_ptr2(chunk_to_used_packets, i, 0),
                                     dims_chunk_to_used_packets_1,
                                     (bool*)fast_array_ptr2(chunk_to_used_packets, j, 0));
            }
        }
    }

    // Backward elimination
    for (npy_intp col = dims_a_1 - 1; col >= 0; col--) {
        if (dirty_rows[col]) continue;

        for (npy_intp row = col - 1; row >= 0; row--) {
            if (dirty_rows[row]) continue;

            if (*((bool*)fast_array_ptr2(A, row, col))) {
                do_xor_bool_optimized((bool*)fast_array_ptr2(A, row, 0),
                                     (bool*)fast_array_ptr2(A, col, 0),
                                     dims_a_1, (bool*)fast_array_ptr2(A, row, 0));
                do_xor_byte_optimized((BYTE*)fast_array_ptr2(b, row, 0),
                                     (BYTE*)fast_array_ptr2(b, col, 0),
                                     dims_b_1, (BYTE*)fast_array_ptr2(b, row, 0));
                do_xor_bool_optimized((bool*)fast_array_ptr2(chunk_to_used_packets, row, 0),
                                     (bool*)fast_array_ptr2(chunk_to_used_packets, col, 0),
                                     dims_chunk_to_used_packets_1,
                                     (bool*)fast_array_ptr2(chunk_to_used_packets, row, 0));
            }
        }
    }

    free(dirty_rows);
    return PyBool_FromLong(!dirty);
}


// Method definitions
static char cdnarules_sat_docs[] =
   "microsatellite(text, lengthToLookFor): Finds the maximum microsatellite with the given length!";
static char cdnarules_lseq_docs[] =
   "longestSequenceOfChar(text, character_to_look_for): Finds the maximum sequence of a char (or any char if given '*')!";
static char cdnarules_rreg_docs[] =
   "repeatRegion(text, lengthToLookFor): returns true if a region of 'lengthToLookFor' is repeated within text (includes overlapping texts)";
static char cdnarules_srreg_docs[] =
   "smallRepeatRegion(text, lengthToLookFor): returns a error-value based on the number repeats the text contains";
static char cdnarules_getq_docs[] =
   "getQUAT(bit1,bit2): returns the DNA base for the given bits";
static char cdnarules_byte2quats_docs[] =
   "byte2QUATS(byte): Converts a given byte to DNA representation";
static char cdnarules_strcsubstr[] =
   "strContainsSub(text,substr): Returns true if the substr is present in text";
static char bitsSet_docs[] =
    "bitsSet(integer): returns the number of bits set int given integer";
static char grayCode_docs[] =
    "grayCode(integer): create a new int for graycode construction";
static char buildGraySequence_docs[] =
    "buildGraySequence(length, b): build up a grey sequence of length <length> with bit b set";
static char bitSet_docs[] =
    "bitSet(X,b): returns if bit b is set in X";
static char xorarray_docs[] =
    "xor_array(X,Y): returns the xor of the two input arrays";
static char elimination_docs[] =
    "elimination(A,b,packet_mapping, chunk_to_used_packets): performs gaussian elimination on A and b. returns true (for now); chunk_to_used_packets MUST be a square matrix >= max(A[rows], A[cols])";
static char elimination_with_first_row_docs[] =
    "elimination_with_first_row(A,b,packet_mapping, chunk_to_used_packets, first_row_idx=-1): performs gaussian elimination on A and b with optional first row index. If first_row_idx is provided and A[first_row_idx, 0] is True, that row will be used as the first pivot. returns true if solved; chunk_to_used_packets MUST be a square matrix >= max(A[rows], A[cols])";
static char gc_content_docs[] =
    "gc_content(text): returns the percentage of GC in the given text";

static PyMethodDef cdnarules_funcs[] = {
   {"bitsSet", bitsSet, METH_VARARGS, bitsSet_docs},
   {"microsatellite", microsatellite, METH_VARARGS, cdnarules_sat_docs},
   {"longestSequenceOfChar", longestSequenceOfChar, METH_VARARGS, cdnarules_lseq_docs},
   {"repeatRegion", repeatRegion, METH_VARARGS, cdnarules_rreg_docs},
   {"smallRepeatRegion", smallRepeatRegion, METH_VARARGS, cdnarules_srreg_docs},
   {"getQUAT", getQUAT, METH_VARARGS, cdnarules_getq_docs},
   {"byte2QUATS", byte2QUATS, METH_VARARGS, cdnarules_byte2quats_docs},
   {"strContainsSub", strContainsSub, METH_VARARGS, cdnarules_strcsubstr},
   {"grayCode", grayCode, METH_VARARGS, grayCode_docs},
   {"buildGraySequence", buildGraySequence, METH_VARARGS, buildGraySequence_docs},
   {"bitSet", bitSet, METH_VARARGS, bitSet_docs},
   {"xorArray", xor_array, METH_VARARGS, xorarray_docs},
   {"elimination", elimination, METH_VARARGS, elimination_docs},
   {"elimination_with_first_row", elimination_with_first_row, METH_VARARGS, elimination_with_first_row_docs},
   {"gc_content", gc_content, METH_VARARGS, gc_content_docs},
   {NULL, NULL, 0, NULL}
};

static struct PyModuleDef cdnarules = {
    PyModuleDef_HEAD_INIT,
    "cdnarules",
    "Extension module for fast DNARules processing with optimizations!",
    -1,
    cdnarules_funcs
};

PyMODINIT_FUNC PyInit_cdnarules(void) {
    import_array();
    return PyModule_Create(&cdnarules);
}