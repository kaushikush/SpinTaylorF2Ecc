# check SWIG Python module wrapping the LAL library
# Author: Karl Wette, 2011

import os
import datetime

def msg(str):
    print(os.path.basename(__file__) + ": " + str)
class error(Exception):
    def __init__(self, str):
        self.str = str
    def __str__(self):
        return str

# check memory allocation
if not cvar.swiglal_debug:
    msg("skipping memory allocation")
else:
    LALCheckMemoryLeaks()
    mem1 = LALDetector()
    mem2 = LALStringVector()
    mem3 = COMPLEX8Vector()
    mem4 = XLALCreateREAL8Vector(3)
    msg("*** below should be an error message from LALCheckMemoryLeaks() ***")
    try:
        LALCheckMemoryLeaks()
        raise error("expected exception")
    except:
        pass
    msg("*** above should be an error message from LALCheckMemoryLeaks() ***")
    del mem1, mem2, mem3
    XLALDestroyREAL8Vector(mem4)
    LALCheckMemoryLeaks()
    msg("passed memory allocation")

# check complex number conversions
assert(XLALCOMPLEX8Add(complex(1, 3), complex(2, -5)) == complex(3, -2))
assert(XLALCOMPLEX8Sub(complex(4, 2), complex(10, 5)) == complex(-6, -3))
assert(XLALCOMPLEX16Mul(complex(10, -7), complex(30, 17)) == complex(419, -40))
assert(XLALCOMPLEX16Div(complex(111.75, -120.50), complex(5, 12)) == complex(-5.25, -11.5))
msg("passed complex number conversions")

# check string conversions
strs = ["a", "bc", "def"]
sv = XLALCreateStringVector(*strs)
assert(sv.length == 3)
assert((sv.data == strs).all())
strs[0] = "ghijk"
sv.data_setel(0, strs[0])
strs.append("lmnopq")
XLALAppendString2Vector(sv, strs[3])
assert(sv.length == 4)
for i in range(0,4):
    assert(sv.data_getel(i) == strs[i])
XLALDestroyStringVector(sv)
msg("passed string conversions")

# check static vector/matrix conversions
if not cvar.swiglal_debug:
    msg("skipping static vector/matrix conversions")
else:
    sts = swiglal_static_test_struct()
    assert(len(sts.vector) == 3)
    assert(len(sts.enum_vector) == 3)
    assert(sts.matrix.shape == (2, 3))
    assert(sts.enum_matrix.shape == (2, 3))
    sts.vector = [3, 2, 1]
    assert((sts.vector == [3, 2, 1]).all())
    sts.matrix = [[4, 5, 6], (9, 8, 7)]
    try:
        sts.matrix = [[1.1, 2.3, 4.5], [6.5, 4.3, 2.1]]
        raise error("expected exception")
    except:
        pass
    assert((sts.matrix == [[4, 5, 6], [9, 8, 7]]).all())
    for i in range(0,3):
        sts.enum_vector_setel(i, 2*i + 3)
        assert(sts.enum_vector_getel(i) == (2*i + 3))
    del sts
    assert(not cvar.swiglal_static_test_vector.any())
    assert(not cvar.swiglal_static_test_matrix.any())
    assert(not cvar.swiglal_static_test_enum_vector.any())
    assert(not cvar.swiglal_static_test_enum_matrix.any())
    cvar.swiglal_static_test_vector = cvar.swiglal_static_test_const_vector
    assert((cvar.swiglal_static_test_vector == [1, 2, 4]).all())
    assert(swiglal_static_test_const_vector_getel(2) == 4)
    try:
        swiglal_static_test_const_vector_getel(20)
        raise error("expected exception")
    except:
        pass
    msg("passed static vector/matrix conversions")

# passed all tests!
msg("================")
msg("PASSED all tests")
msg("================")
