# Define Task object
module Psi
  class Task

    def run
      # IO.read from Ruby standard library
      # ruby_file define by Psi4 to be a member of Task
      # see task.cc:create_ruby_class(...)
      file = IO.read(ruby_file)
      eval file
    end
    
  end
end
