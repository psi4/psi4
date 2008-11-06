# Define Task object
module Psi
  class Task

    def run
      # IO.read from Ruby standard library
      # input define by Psi4 to be a member of Task
      # see task.cc:create_ruby_class(...)
      file = IO.read(input)
      eval file
    end
    
  end
end
