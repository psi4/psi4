import pytest
import psi4

pytestmark = [pytest.mark.psi, pytest.mark.api]

def test_option_group():
    options = psi4.core.FOptions()
    groupname = "groupname"
    options.set_group(groupname)
    assert options.get_group() == groupname

def test_option_addgetset():
    options = psi4.core.FOptions()

    # None
    options.add_bool("key1", None)
    options.add_bool("key2", None, "description")
    options.set_bool("key1", True)
    assert options.get_bool("key1") is True

    # bool
    options.add_bool("key3", True)
    options.add_bool("key4", False, "description")
    options.add_bool("key5", 42)
    options.add_bool("key6", [], "description")
    assert options.get_bool("key3") is True
    assert options.get_bool("key4") is False
    assert options.get_bool("key5") is True
    assert options.get_bool("key6") is False
    
    # int
    options.add_int("key_int", 42)
    assert options.get_int("key_int") == 42
    options.set_int("key_int", 43)
    assert options.get_int("key_int") == 43
    
    # double
    options.add_double("key_double", 42.0)
    assert options.get_double("key_double") == 42.0
    options.set_double("key_double", 43.0)
    assert options.get_double("key_double") == 43.0

    # str
    options.add_str("key_str", "foo")
    assert options.get_str("key_str") == "FOO"
    options.set_str("key_str", "bar")
    assert options.get_str("key_str") == "BAR"
    
    options.add_str("key_str2", None, ["DCT", "OCC", "SCF"])
    options.set_str("key_str2", "dct")
    assert options.get_str("key_str2") == "DCT"

    # int list
    options.add_int_list("key_int_list")
    options.set_int_list("key_int_list", [1, 2])
    assert options.get_int_list("key_int_list") == [1, 2]

    # double list
    options.add_double_list("key_double_list")
    options.set_double_list("key_double_list", [1.0, 2.0])
    assert options.get_double_list("key_double_list") == [1.0, 2.0]

def test_option_errors():
    options = psi4.core.FOptions()

    # can't get/set unregistered keys
    with pytest.raises(RuntimeError) as err:
        options.get_bool("key1")
    
    assert 'not registered' in str(err.value)
    with pytest.raises(RuntimeError) as err:
        options.set_bool("key1", True)
    
    assert 'not registered' in str(err.value)

    options.add_int("key_int", None)

    # can't duplicate keys
    options.add_bool("key1", False)
    with pytest.raises(RuntimeError) as err:
        options.add_bool("key1", False)
    
    assert 'already used' in str(err.value)

    with pytest.raises(RuntimeError) as err:
        options.add_bool("key1", False, "description")
    
    assert 'already used' in str(err.value)
    with pytest.raises(RuntimeError) as err:
        options.add_int("key1", 42)
    
    assert 'already used' in str(err.value)

    # can't get none
    with pytest.raises(RuntimeError) as err:
        options.get_int("key_int")
    
    assert 'is set to None' in str(err.value)

    # getters must be of correct type
    with pytest.raises(RuntimeError) as err:
        options.get_int("key1")
    
    assert 'type for this option is' in str(err.value)
    with pytest.raises(RuntimeError) as err:
        options.get_double("key1")
    
    assert 'type for this option is' in str(err.value)
    with pytest.raises(RuntimeError) as err:
        options.get_str("key1")
    
    assert 'type for this option is' in str(err.value)

    # setters must be of correct type
    with pytest.raises(RuntimeError) as err:
        options.set_int("key1", 42)
    
    assert 'type for this option is' in str(err.value)
    with pytest.raises(RuntimeError) as err:
        options.set_double("key1", 42.0)
    
    assert 'type for this option is' in str(err.value)
    with pytest.raises(RuntimeError) as err:
        options.set_str("key1", "foo")
    
    assert 'type for this option is' in str(err.value)

    # str must be in allowed values
    options.add_str("key_str", None, ["DCT", "OCC", "SCF"])
    with pytest.raises(RuntimeError) as err:
        options.set_str("key_str", "dfocc")
    
    assert 'in allowed_values' in str(err.value)

