import unittest

import sys
import os

# Add the directory containing the 'grendel_tests' package to sys.path
from grendel.util.containers import AliasedDict
from grendel_tests import expected_failure

sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))

class AliasedDictTest(unittest.TestCase):

    def setUp(self):
        self.d = AliasedDict(
            {
                ("hello", "hello world", 17) : 170,
                ("goodbye", "goodbye world", 75, 31) : "the end",
                frozenset({"failure", "is not", "success"}) : 3.14,
                8 : 88,
                frozenset((9, 10)): [1, 2, 3],
                170 : 17
            }
        )

    def tearDown(self):
        del self.d

    #region | Constructor errors |

    def test_construct_error_1(self):
        with self.assertRaises(KeyError):
            AliasedDict({ ("a", "b") : 3, "a" : 5})

    def test_construct_error_2(self):
        with self.assertRaises(KeyError):
            AliasedDict({ ("a", "b") : 3, ("b", "c") : 4})

    def test_construct_error_3(self):
        with self.assertRaises(KeyError):
            AliasedDict({ ("a", ("b", "c")) : 3})

    def test_construct_error_4(self):
        with self.assertRaises(KeyError):
            AliasedDict({ ("a", frozenset(("b", "c"))) : 3})

    #endregion

    def test_len(self):
        self.assertEqual(len(self.d), 6)

    #region | get |

    def test_get1(self):
        self.assertEqual(self.d['hello'], 170)

    def test_get2(self):
        self.assertEqual(self.d[17], 170)

    def test_get3(self):
        self.assertEqual(self.d[8], 88)

    def test_get4(self):
        self.assertEqual(self.d[10], [1, 2, 3])

    def test_get5(self):
        self.assertEqual(self.d[170], 17)

    def test_get6(self):
        self.assertEqual(self.d[{"failure", "is not", "success"}], 3.14)

    def test_get7(self):
        self.assertEqual(self.d[("is not", "failure", "success")], 3.14)

    def test_get8(self):
        self.assertEqual(self.d[frozenset(("is not", "failure", "success"))], 3.14)

    def test_get9(self):
        self.d[10][1] = 13
        self.assertEqual(self.d[9], [1, 13, 3])

    def test_get_missing_1(self):
        with self.assertRaises(KeyError):
            self.d[32]

    def test_get_missing_2(self):
        with self.assertRaises(KeyError):
            self.d[("hello", 17)]

    def test_get_missing_3(self):
        with self.assertRaises(KeyError):
            self.d[(8, 9)]

    #endregion

    #region | set |

    def test_set1(self):
        self.d[1] = 2
        self.assertEqual(self.d[1], 2)
        self.assertEqual(len(self.d), 7)

    def test_set2(self):
        self.d["hello"] = 13
        self.assertEqual(self.d[17], 13)
        self.assertEqual(len(self.d), 6)

    def test_set3(self):
        self.d[(9, 10)] = 7, 8
        self.assertEqual(self.d[9], (7, 8))
        self.assertEqual(len(self.d), 6)

    def test_set4(self):
        self.d[(13, 14)] = 15
        self.assertEqual(self.d[14], 15)
        self.assertEqual(len(self.d), 7)

    def test_set5(self):
        self.d[9, 10] = 7, 8
        self.assertEqual(self.d[9], (7, 8))
        self.assertEqual(len(self.d), 6)

    def test_set_error_1(self):
        with self.assertRaises(KeyError):
            self.d[("failure", "is not")] = 5

    def test_set_error_2(self):
        with self.assertRaises(KeyError):
            self.d["failure", 3] = 5

    def test_set_error_3(self):
        with self.assertRaises(KeyError):
            self.d[170, 171] = 5

    def test_set_error_4(self):
        with self.assertRaises(KeyError):
            self.d[(170, 17)] = 5

    def test_set_error_5(self):
        with self.assertRaises(TypeError):
            self.d[list((1,2,3,4))] = 5

    def test_set_error_6(self):
        with self.assertRaises(KeyError):
            self.d[1,(2,3),4] = 5

    #endregion

    #region | del |

    def test_del1(self):
        del self.d[170]
        self.assertEqual(len(self.d), 5)

    def test_del2(self):
        del self.d["failure"]
        self.assertEqual(len(self.d), 5)

    def test_del3(self):
        del self.d["failure", "is not", "success"]
        self.assertEqual(len(self.d), 5)

    def test_del_error_1(self):
        with self.assertRaises(KeyError):
            del self.d["failure", "is not"]

    def test_del_error_2(self):
        with self.assertRaises(KeyError):
            del self.d[1]

    #endregion

    #region | contains |
    def test_contains_1(self):
        self.assertIn("hello", self.d)

    def test_contains_2(self):
        self.assertIn({"hello", "hello world", 17}, self.d)

    def test_contains_3(self):
        self.assertIn(8, self.d)

    def test_contains_4(self):
        self.assertIn((8,), self.d)

    def test_contains_5(self):
        self.assertNotIn(171, self.d)

    def test_contains_6(self):
        self.assertNotIn(("goodbye", "goodbye world", 75), self.d)

    def test_contains_error(self):
        with self.assertRaises(TypeError):
            self.assertNotIn([9, 10], self.d)

    #endregion

    # TODO More unit tests

    #region | update |

    def test_update_1(self):
        self.d.update({10 : 12})
        self.assertEqual(self.d[9], 12)

    def test_update_2(self):
        self.d.update(hello=3)
        self.assertEqual(self.d["hello world"], 3)

    def test_update_3(self):
        # Matches standard library behavior
        self.d.update({"hello" : 2}, hello=3)
        self.assertEqual(self.d["hello world"], 3)

    #endregion

    #region | fromkeys |

    def test_fromkeys_1(self):
        d2 = AliasedDict.fromkeys((("hello", "world"), 5, (13, 17)))
        d2["world"] = 3.14
        self.assertEqual(d2["hello"], 3.14)

    def test_fromkeys_2(self):
        d2 = AliasedDict.fromkeys((("hello", "world"), 5, (13, 17)))
        d2[17] = 2
        self.assertEqual(d2[13], 2)

    def test_fromkeys_3(self):
        d2 = AliasedDict.fromkeys((("hello", "world"), 5, (13, 17)))
        d2[5] = 8
        self.assertEqual(d2[5], 8)

    #endregion

    @expected_failure
    def test_viewkeys(self):
        v = self.d.viewkeys()
        self.assertEqual(len(v), 6)

    # End AliasedDictTest
    pass

