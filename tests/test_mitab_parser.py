import unittest

from BioTools.MITAB_parser import Check_Value


class TestCheckValue(unittest.TestCase):
    def test_accepts_valid_value(self):
        Check_Value("human", {"human", "mouse"}, "species")

    def test_raises_for_invalid_value(self):
        with self.assertRaises(Exception):
            Check_Value("yeast", {"human", "mouse"}, "species")


if __name__ == "__main__":
    unittest.main()
