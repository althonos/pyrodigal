from libc.stdint cimport uint32_t


cdef extern from *:

    ctypedef uint32_t Py_UCS4

    object PyUnicode_New(Py_ssize_t, Py_UCS4)
    int PyUnicode_WriteChar(object, Py_ssize_t, Py_UCS4) except -1
    Py_UCS4 PyUnicode_ReadChar(object, Py_ssize_t)

    cdef int PyUnicode_1BYTE_KIND
    cdef int PyUnicode_2BYTE_KIND
    cdef int PyUnicode_4BYTE_KIND

    int     PyUnicode_KIND(object)
    void*   PyUnicode_DATA(object)
    void    PyUnicode_WRITE(int kind, void *data, Py_ssize_t index, Py_UCS4 value) nogil
    Py_UCS4 PyUnicode_READ (int kind, void *data, Py_ssize_t index) nogil

    Py_ssize_t PyUnicode_GET_LENGTH(object)

    object PyUnicode_FromStringAndSize(const char*, Py_ssize_t)
    object PyUnicode_FromKindAndData(int kind, const void *buffer, Py_ssize_t size)

    int    PyUnicode_READY(object) except -1
