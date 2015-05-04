/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */
#ifndef PSIAPI_LABELS_H_
#define PSIAPI_LABELS_H_

namespace PsiAPI{
/** \brief The recognized list of labels
 *
 *  In order to have a "typesafe" way of checking for labels I have
 *  elected for an enumeration over strings.  With strings we would
 *  have to worry about different spellings, capitalization, etc.
 *  which may lead to subtle hard to find errors.  With enums the compiler
 *  will yell at you if you select one that doesn't exist.
 *
 *  For details about the label system see the %Atom class.
 *
 *  This is the comprehensive list of labels and what they mean.  If you
 *  add a label make sure you document it (the appropriate place would be
 *  here fyi).
 *
 *  Label (Description):
 *     - NOLABEL The default label.  Indicates that there is nothing
 *               special about the atom.
 *     - GHOST   The atom is to be treated as a ghost atom (basis functions
 *               with no electrons)
 *     - DUMMY   The atom is used for specify a point in space.
 *     - MM      The atom is part of the MM part of a QM/MM job. Possessing
 *               both an MM and a QM tag indicates belonging to the interface.
 *     - QM      The atom is part of the QM part of the QM/MM job. Possessing
 *               both an MM and a QM tag indicates belonging to the interface.
 *     - CHARGE  The atom is a point charge.
 *
 *
 */
enum AtomLabels{NOLABEL,GHOST,DUMMY,MM,QM,CHARGE};
}



#endif /* PSIAPI_LABELS_H_ */
