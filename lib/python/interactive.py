from psi4 import *
import readline # optional, will allow Up/Down/History in the console
import code

def run():
    print_out("\nStarting interactive session.\n\n")

    vars = globals().copy()
    vars.update(locals())
    shell = code.InteractiveConsole(vars)
    shell.interact()

