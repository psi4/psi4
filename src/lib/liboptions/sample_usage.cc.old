  try {
    Options options;
    
    std::cout << "testing Options class" << std::endl;
    options.add("mult", 1);
    options.add("array", new ArrayType());
    options.add("map", new MapType());
    options.add("reference", "rhf", "rhf rohf uhf tcscf");
    
    options["array"].add(3);                       // 0
    options["array"].add(0);                       // 1
    options["array"].add("3rd", "");               // 2
    options["array"].add(new ArrayType());         // 3
    options["array"][3].add("inner", "");          // 3-0
    options["array"][3].add(new ArrayType());      // 3-1
    options["array"][3][1].add("inner-inner", ""); // 3-1-0
    
    options["map"].add("test", 10);
    options["map"].add("array", new ArrayType());
    options["map"]["array"].add(20);
    
    std::cout << "mult " << options["mult"].to_string() << std::endl;
    std::cout << "array " << options["array"].to_string() << std::endl;
    
    std::cout << "array[2] " << options["array"][2].to_string() << std::endl;
    std::cout << "array[3] " << options["array"][3].to_string() << std::endl;
    
    std::cout << "map " << options["map"].to_string() << std::endl << std::endl;

    // Testing assign
    options["mult"].assign(2);
    options["map"]["array"][0].assign(256);
    
    options.print();
    
    // Test StringDataType choices
    options["reference"].assign("ROHF");
    
    options.print();
  }
  catch (IndexException e) {
    std::cout << "IndexException caught: " << e.what() << std::endl;
  }
  catch (DuplicateKeyException e) {
    std::cout << "DuplicateKeyException caught: " << e.what() << std::endl;
  }
  catch (DataTypeException e) {
    std::cout << "DataTypeException caught: " << e.what() << std::endl;
  }
  catch (NotImplementedException e) {
    std::cout << "NotImplementedException caught: " << e.what() << std::endl;
  }

