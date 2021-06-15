from libc.stdint cimport uint32_t


cdef extern from *:

    ctypedef uint32_t Py_UCS4

    object PyUnicode_New(Py_ssize_t, Py_UCS4)
    int PyUnicode_WriteChar(object, Py_ssize_t, Py_UCS4) except -1
    Py_UCS4 PyUnicode_ReadChar(object, Py_ssize_t)
