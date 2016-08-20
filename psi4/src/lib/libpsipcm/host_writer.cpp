/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "psi4/psi4-dec.h" //Gives us psi::outfile

extern "C" void host_writer(const char * message, int /* message_length */)
{
  psi::outfile->Printf(message);
  psi::outfile->Printf("\n");
}
