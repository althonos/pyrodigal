from libc.stdint cimport uint8_t, uint16_t, uint32_t


cdef extern from *:

    ctypedef uint8_t  Py_UCS1
    ctypedef uint16_t Py_UCS2
    ctypedef uint32_t Py_UCS4

    int PyUnicode_KIND(object)
    cdef int PyUnicode_1BYTE_KIND
    cdef int PyUnicode_2BYTE_KIND
    cdef int PyUnicode_4BYTE_KIND

    Py_UCS1* PyUnicode_1BYTE_DATA(object)
    Py_UCS2* PyUnicode_2BYTE_DATA(object)
    Py_UCS4* PyUnicode_4BYTE_DATA(object)

    object PyUnicode_New(Py_ssize_t, Py_UCS4)
    int PyUnicode_WriteChar(object, Py_ssize_t, Py_UCS4)
