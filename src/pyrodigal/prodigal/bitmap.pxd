# a typedef to make it more obvious when a sequence is a string (char*) and
# when it is a bitmap (unsigned char*)
ctypedef unsigned char* bitmap_t


cdef extern from "bitmap.h" nogil:

    unsigned char test(bitmap_t, int) noexcept
    void clear(bitmap_t, int) noexcept
    void set(bitmap_t, int) noexcept
    void toggle(bitmap_t, int) noexcept
