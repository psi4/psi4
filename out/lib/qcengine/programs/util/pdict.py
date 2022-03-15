from decimal import ROUND_CEILING, ROUND_FLOOR, Decimal

import numpy as np

# This file has considerable history in the quantum chemistry package Psi4


class PreservingDict(dict):
    """Class to store quantum chemical quantities extracted from output
    files. Extends the dictionary object to (1) store key as all-caps
    version of itself and (2) validate value for duplicate values for the
    same key by testing which has more decimal places and whether value
    the same within a plausible rounding error. Allows consistency checks
    when parsing output files without loss of precision.

    works in decimal.Decimal (scalar) and np.ndarray (non-scalar)

    """

    def __init__(self, *args, **kwargs):
        self.update(*args, **kwargs)

    def __setitem__(self, key, value, accept_places=11):
        try:
            key = key.upper()
        except AttributeError:
            raise AttributeError("Keys stored as upper-case strings: %s unsuitable" % (key))

        if isinstance(value, (list, np.ndarray)):
            # non-scalar
            value = np.array(value)

            if key in self.keys() and not "CURRENT" in key:
                if np.allclose(self[key], value, atol=1.0e-5, rtol=1e-3):
                    best_value = value  # no way to choose, really
                else:
                    raise ValueError(
                        """Output file yielded both {} and {} as values for quantity {}.""".format(
                            self[key], value, key
                        )
                    )
                # print("""Resetting array {} to {}""".format(key, best_value))
            else:
                best_value = value
                # print("""Setting   array {} to {}""".format(key, best_value))

        else:
            # scalar
            value = Decimal(value)
            if abs(float(value)) < 1.0e-16:
                value = Decimal("0")

            if key in self.keys() and not "CURRENT" in key:
                # Validate choosing more detailed value for variable
                existing_exp = self[key].as_tuple().exponent  # 0.1111 --> -4
                candidate_exp = value.as_tuple().exponent
                if existing_exp > candidate_exp:  # candidate has more digits
                    places = Decimal(10) ** (existing_exp + 1)  # exp+1 permits slack in rounding
                    best_value = value
                else:  # existing has more digits
                    places = Decimal(10) ** (candidate_exp + 1)
                    best_value = self[key]
                # DEBUG print(f"{existing_exp=} {candidate_exp=} {places=} {self[key]=} {value=}")
                # Validate values are the same
                # places = max(places, Decimal('1E-11'))  # for computed psivars
                places = max(places, Decimal("1E-{}".format(accept_places)))  # for computed psivars
                # print('FLOOR: ', self[key].quantize(places, rounding=ROUND_FLOOR) - value.quantize(places, rounding=ROUND_FLOOR))
                # print('CEIL:  ', self[key].quantize(places, rounding=ROUND_CEILING) - value.quantize(places, rounding=ROUND_CEILING))
                if (
                    self[key]
                    .quantize(places, rounding=ROUND_CEILING)
                    .compare(value.quantize(places, rounding=ROUND_CEILING))
                    != 0
                ) and (
                    self[key]
                    .quantize(places, rounding=ROUND_FLOOR)
                    .compare(value.quantize(places, rounding=ROUND_FLOOR))
                    != 0
                ):
                    raise ValueError(
                        """Output file yielded both %s and %s as values for quantity %s."""
                        % (self[key].to_eng_string(), value.to_eng_string(), key)
                    )
                # print("""Resetting variable {} to {}""".format(key, best_value.to_eng_string()))
            else:
                best_value = value
                # print("""Setting   variable {} to {}""".format(key, best_value.to_eng_string()))

        super(PreservingDict, self).__setitem__(key, best_value)

    def update(self, *args, **kwargs):
        if args:
            if len(args) > 1:
                raise TypeError("update expected at most 1 arguments, " "got %d" % len(args))
            other = dict(args[0])
            for key in other:
                self[key] = other[key]
        for key in kwargs:
            self[key] = kwargs[key]

    def setdefault(self, key, value=None):
        if key not in self:
            self[key] = value
        return self[key]


if __name__ == "__main__":
    c4info = PreservingDict()
    c4info["scf 4.5e0 total energy"] = "-1.e-4"
    c4info["1.3"] = ".4"
    c4info["curl"] = "-437.12345678"
    c4info["curl"] = "-437.12345677"
    c4info["curl"] = "-437.123456"
    c4info["curl"] = "-437.123457"
    c4info["curl"] = "-437.1234444"  # fails
    c4info["curl"] = "-437.123456789"
    c4info["curl"] = "-437.1234567779"  # fails
    print(c4info)
