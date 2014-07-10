/** 
    Creation, cloning, and destruction of array distribution templates and 
    distributed array descriptors through two CCA ports.  

    Templates conform to the DistArrayTemplate abstract base class,
    and descriptors to the DistArrayDescriptor abstract base class.  

    In this implementation, some care is taken to try to check for bad
    pointers and other problems in order to return an error rather
    than cause the code to seg fault.
*/

#include <string>
#include <stdlib.h>
#include <stdio.h>

// Our header
#include "gacca.h"

// Concrete classes for templates, descriptors and distributed arrays
#include "DADFTemplate.h"
#include "DADFDescriptor.h"
#include "GA_DADFArray.h"

/******************************************************************************
 * DistArrayTemplFactoryPort
******************************************************************************/

/** Return an uninitialized template object. */
DistArrayTemplate*
GA::GAServices::createTemplate(std::string name) {
  DADFTemplate* templDADF;

  /** Create a template of our private class */

  templDADF = new DADFTemplate( name );
  if ( templDADF == 0 ) {
    cerr << "GA::GAServices::createTemplate:" <<
      "Unable to create new template." << endl;
    return 0;
  }

  /** Record the address for later reference */

  templateRecord[ templDADF ]++;

  /** Cast from our private class to the abstract base class we
      publicized */

  return dynamic_cast<DistArrayTemplate *> (templDADF);
}

/** Return a template object initialized with the contents of
    another, but not frozen against modification. */
DistArrayTemplate* 
GA::GAServices::cloneTemplate(
  DistArrayTemplate * original, std::string cloneName) {

  DADFTemplate * originalDADF, * templDADF;

  /** Check that the original template is one we provided.  */

  if ( templateRecord.count( original ) ) {

    /** How many templates do we have with this address?
	0 means we had one, but its been destroyed already. */

    if ( templateRecord[ original ] == 0 ) {
      cerr << "GA::GAServices::cloneTemplate: " <<
	"Original template invalid: previously destroyed." << endl;
      return 0;
    } 

    /** Anything other than 1 shouldn't happen.  We'll let it go with
	a warning -- in case it is our mistake -- but it will probably
	bomb. Maybe this should be an error? */
    else if ( templateRecord[ original ] != 1 ) {
      cerr << "GA::GAServices::cloneTemplate: " <<
	"Internal records indicate " << templateRecord[ original ] <<
	"outstanding templates at address of original." << endl <<
	"Something has gone badly wrong. Will try to proceed anyway."
	   << endl;
    }
  } 

  /** If we don't even have a record of this template, it means we
      didn't create it.  So either its a bad pointer, or it was
      created with another factory.  In either case, we don't really
      want to mess with it. */

  else {
    cerr << "GA::GAServices::cloneTemplate: " <<
      "Original template not from this factory." << endl;
    return 0;
  }

  /** Cast it from the abstract class to our private class.  This is a
      further check that we should be able to deal with it properly,
      and something we have to do to create a new instance. */

  originalDADF = dynamic_cast<DADFTemplate *>(original);
  if ( originalDADF == 0 ) {
    cerr << "GA::GAServices::cloneTemplate: " <<
      "Original template is of the wrong type." << endl;
    return 0;
  }

  /** Create a new instance with a copy constructor from the original */

  templDADF = new DADFTemplate( cloneName, *originalDADF );
  if ( templDADF == 0 ) {
    cerr << "GA::GAServices::createTemplate: " <<
      "Unable to create new template." << endl;
    return 0;
  }

  /** Record the address for our records */

  templateRecord[ templDADF ]++;

  /** Cast from our private class to the abstract base class we
      publicized */

  return dynamic_cast<DistArrayTemplate *> (templDADF);
}

/** Destroy an existing template. 

    @returns a null pointer in victim 
*/
int GA::GAServices::destroyTemplate(
  DistArrayTemplate * & victim) { 
  DADFTemplate * victimDADF;

  /** Check that the victim template is one we provided.  */

  if ( templateRecord.count( victim ) ) {

    /** How many templates do we have with this address?
	0 means we had one, but its been destroyed already. */

    if ( templateRecord[ victim ] == 0 ) {
      cerr << "GA::GAServices::destroyTemplate: " <<
	"Template invalid: previously destroyed." << endl;
      return -2;
    } 

    /** Anything other than 1 shouldn't happen.  We'll let it go with
	a warning -- in case it is our mistake -- but it will probably
	bomb. Maybe this should be an error? */
    else if ( templateRecord[ victim ] != 1 ) {
      cerr << "GA::GAServices::destroyTemplate: " <<
	"Internal records indicate " << templateRecord[ victim ] <<
	"outstanding templates at address of victim." << endl <<
	"Something has gone badly wrong. Will try to proceed anyway."
	   << endl;
    }
  } 

  /** If we don't even have a record of this template, it means we
      didn't create it.  So either its a bad pointer, or it was
      created with another factory.  In either case, we don't really
      want to mess with it. */

  else {
    cerr << "GA::GAServices::destroyTemplate: " <<
      "Template not from this factory." << endl;
    return -1;
  }

  /** Cast it from the abstract class to our private class.  This is a
      further check that we should be able to deal with it properly. */

  victimDADF = dynamic_cast<DADFTemplate *>(victim);
  if ( victimDADF == 0 ) {
    cerr << "GA::GAServices::destroyTemplate: " <<
      "Template is of the wrong type." << endl;
    return -3;
  }

  /** Get rid of the victim */
  delete victimDADF;

  /** Record the fact that we got rid of it */
  templateRecord[ victim ]--;

  /** Clear the pointer so user knows its no longer valid */
  victim = 0;

  /** Success */
  return 0;
}

/******************************************************************************
 * DistArrayDescrFactoryPort
******************************************************************************/

/** Return an uninitialized descriptor object. */
DistArrayDescriptor*
GA::GAServices::createDescriptor(std::string name) {
  DADFDescriptor* descrDADF;

  /** Create a descriptor of our private class */

  descrDADF = new DADFDescriptor( name );
  if ( descrDADF == 0 ) {
    cerr << "GA::GAServices::createDescriptor:" <<
      "Unable to create new descriptor." << endl;
    return 0;
  }

  /** Record the address for later reference */

  descriptorRecord[ descrDADF ]++;

  /** Cast from our private class to the abstract base class we
      publicized */

  return dynamic_cast<DistArrayDescriptor *> (descrDADF);
}

/** Return a descriptor object initialized with the contents of
    another, but not frozen against modification. */
DistArrayDescriptor* 
GA::GAServices::cloneDescriptor(
  DistArrayDescriptor * original, std::string cloneName) {

  DADFDescriptor * originalDADF, * descrDADF;

  /** Check that the original descriptor is one we provided.  */

  if ( descriptorRecord.count( original ) ) {

    /** How many descriptors do we have with this address?
	0 means we had one, but its been destroyed already. */

    if ( descriptorRecord[ original ] == 0 ) {
      cerr << "GA::GAServices::cloneDescriptor: " <<
	"Original descriptor invalid: previously destroyed." << endl;
      return 0;
    } 

    /** Anything other than 1 shouldn't happen.  We'll let it go with
	a warning -- in case it is our mistake -- but it will probably
	bomb. Maybe this should be an error? */
    else if ( descriptorRecord[ original ] != 1 ) {
      cerr << "GA::GAServices::cloneDescriptor: " <<
	"Internal records indicate " << descriptorRecord[ original ] <<
	"outstanding descriptors at address of original." << endl <<
	"Something has gone badly wrong. Will try to proceed anyway."
	   << endl;
    }
  } 

  /** If we don't even have a record of this descriptor, it means we
      didn't create it.  So either its a bad pointer, or it was
      created with another factory.  In either case, we don't really
      want to mess with it. */

  else {
    cerr << "GA::GAServices::cloneDescriptor: " <<
      "Original descriptor not from this factory." << endl;
    return 0;
  }

  /** Cast it from the abstract class to our private class.  This is a
      further check that we should be able to deal with it properly,
      and something we have to do to create a new instance. */

  originalDADF = dynamic_cast<DADFDescriptor *>(original);
  if ( originalDADF == 0 ) {
    cerr << "GA::GAServices::cloneDescriptor: " <<
      "Original descriptor is of the wrong type." << endl;
    return 0;
  }

  /** Create a new instance with a copy constructor from the original */

  descrDADF = new DADFDescriptor( cloneName, *originalDADF );
  if ( descrDADF == 0 ) {
    cerr << "GA::GAServices::createDescriptor: " <<
      "Unable to create new descriptor." << endl;
    return 0;
  }

  /** Record the address for our records */

  descriptorRecord[ descrDADF ]++;

  /** Cast from our private class to the abstract base class we
      publicized */

  return dynamic_cast<DistArrayDescriptor *> (descrDADF);
}

/** Destroy an existing descriptor. 

    @returns a null pointer in victim
*/
int GA::GAServices::destroyDescriptor(
  DistArrayDescriptor * & victim) { 
  DADFDescriptor * victimDADF;

  /** Check that the victim descriptor is one we provided.  */

  if ( descriptorRecord.count( victim ) ) {

    /** How many descriptors do we have with this address?
	0 means we had one, but its been destroyed already. */

    if ( descriptorRecord[ victim ] == 0 ) {
      cerr << "GA::GAServices::destroyDescriptor: " <<
	"Descriptor invalid: previously destroyed." << endl;
      return -2;
    } 

    /** Anything other than 1 shouldn't happen.  We'll let it go with
	a warning -- in case it is our mistake -- but it will probably
	bomb. Maybe this should be an error? */
    else if ( descriptorRecord[ victim ] != 1 ) {
      cerr << "GA::GAServices::destroyDescriptor: " <<
	"Internal records indicate " << descriptorRecord[ victim ] <<
	"outstanding descriptors at address of victim." << endl <<
	"Something has gone badly wrong. Will try to proceed anyway."
	   << endl;
    }
  } 

  /** If we don't even have a record of this descriptor, it means we
      didn't create it.  So either its a bad pointer, or it was
      created with another factory.  In either case, we don't really
      want to mess with it. */

  else {
    cerr << "GA::GAServices::destroyDescriptor: " <<
      "Descriptor not from this factory." << endl;
    return -1;
  }

  /** Cast it from the abstract class to our private class.  This is a
      further check that we should be able to deal with it properly. */

  victimDADF = dynamic_cast<DADFDescriptor *>(victim);
  if ( victimDADF == 0 ) {
    cerr << "GA::GAServices::destroyDescriptor: " <<
      "Descriptor is of the wrong type." << endl;
    return -3;
  }

  /** Get rid of the victim */
  delete victimDADF;

  /** Record the fact that we got rid of it */
  descriptorRecord[ victim ]--;

  /** Clear the pointer so user knows its no longer valid */
  victim = 0;

  /** Success */
  return 0;
}


/** Return an uninitialized distributed array object. */
DistArray*
GA::GAServices::createArray(std::string name) {
  DADFArray* arrayDADF;

  /** Create an array of our private class */

  arrayDADF = new DADFArray( name );
  if ( arrayDADF == 0 ) {
    cerr << "GA::GAServices::createArray:" <<
      "Unable to create new array." << endl;
    return 0;
  }

  /** Record the address for later reference */
  arrayRecord[ arrayDADF ]++;

  /** Cast from our private class to the abstract base class we
      publicized */

  return dynamic_cast<DistArray *> (arrayDADF);
}

/** Return an array object initialized with the contents of
    another, but not frozen against modification. */
DistArray* 
GA::GAServices::cloneArray(
  DistArray * original, std::string cloneName) {

  DADFArray * originalDADF, * arrayDADF;

  /** Check that the original array is one we provided.  */

  if ( arrayRecord.count( original ) ) {

    /** How many arrays do we have with this address?
	0 means we had one, but its been destroyed already. */

    if ( arrayRecord[ original ] == 0 ) {
      cerr << "GA::GAServices::cloneArray: " <<
	"Original array invalid: previously destroyed." << endl;
      return 0;
    } 

    /** Anything other than 1 shouldn't happen.  We'll let it go with
	a warning -- in case it is our mistake -- but it will probably
	bomb. Maybe this should be an error? */
    else if ( arrayRecord[ original ] != 1 ) {
      cerr << "GA::GAServices::cloneArray: " <<
	"Internal records indicate " << arrayRecord[ original ] <<
	"outstanding arrays at address of original." << endl <<
	"Something has gone badly wrong. Will try to proceed anyway."
	   << endl;
    }
  } 

  /** If we don't even have a record of this array, it means we
      didn't create it.  So either its a bad pointer, or it was
      created with another factory.  In either case, we don't really
      want to mess with it. */

  else {
    cerr << "GA::GAServices::cloneArray: " <<
      "Original array not from this factory." << endl;
    return 0;
  }

  /** Cast it from the abstract class to our private class.  This is a
      further check that we should be able to deal with it properly,
      and something we have to do to create a new instance. */

  originalDADF = dynamic_cast<DADFArray *>(original);
  if ( originalDADF == 0 ) {
    cerr << "GA::GAServices::cloneArray: " <<
      "Original array is of the wrong type." << endl;
    return 0;
  }

  /** Create a new instance with a copy constructor from the original */

  arrayDADF = new DADFArray( cloneName, *originalDADF );
  if ( arrayDADF == 0 ) {
    cerr << "GA::GAServices::createArray: " <<
      "Unable to create new array." << endl;
    return 0;
  }

  /** Record the address for our records */

  arrayRecord[ arrayDADF ]++;

  /** Cast from our private class to the abstract base class we
      publicized */

  return dynamic_cast<DistArray *> (arrayDADF);
}

/** Destroy an existing array. 

    @returns a null pointer in victim
*/
int GA::GAServices::destroyArray(
  DistArray * & victim) { 
  DADFArray * victimDADF;

  /** Check that the victim array is one we provided.  */

  if ( arrayRecord.count( victim ) ) {

    /** How many arrays do we have with this address?
	0 means we had one, but its been destroyed already. */

    if ( arrayRecord[ victim ] == 0 ) {
      cerr << "GA::GAServices::destroyArray: " <<
	"Array invalid: previously destroyed." << endl;
      return -2;
    } 

    /** Anything other than 1 shouldn't happen.  We'll let it go with
	a warning -- in case it is our mistake -- but it will probably
	bomb. Maybe this should be an error? */
    else if ( arrayRecord[ victim ] != 1 ) {
      cerr << "GA::GAServices::destroyArray: " <<
	"Internal records indicate " << arrayRecord[ victim ] <<
	"outstanding arrays at address of victim." << endl <<
	"Something has gone badly wrong. Will try to proceed anyway."
	   << endl;
    }
  } 

  /** If we don't even have a record of this array, it means we
      didn't create it.  So either its a bad pointer, or it was
      created with another factory.  In either case, we don't really
      want to mess with it. */

  else {
    cerr << "GA::GAServices::destroyArray: " <<
      "Array not from this factory." << endl;
    return -1;
  }

  /** Cast it from the abstract class to our private class.  This is a
      further check that we should be able to deal with it properly. */

  victimDADF = dynamic_cast<DADFArray *>(victim);
  if ( victimDADF == 0 ) {
    cerr << "GA::GAServices::destroyArray: " <<
      "Array is of the wrong type." << endl;
    return -3;
  }

  /** Get rid of the victim */
  delete victimDADF;

  /** Record the fact that we got rid of it */
  arrayRecord[ victim ]--;

  /** Clear the pointer so user knows its no longer valid */
  victim = 0;

  /** Success */
  return 0;
}

/******************************************************************************
 * Private stuff
******************************************************************************/

/** List outstanding templates on stderr

    @param label (in) String to label output 
*/
void GA::GAServices::listTemplates(std::string label) {
  std::map<DistArrayTemplate *, int>::iterator it = templateRecord.begin();

  for ( ; it != templateRecord.end(); ++it) {
    if ( it->second == 1) {
      cerr << label << ": " <<
	"WARNING! template `" << (it->first)->getName() <<
	"' at address 0x" << it->first << " still exists." << endl;
    }

    /** Report any bizarrities */

    else if ( it->second < 0 || it->second > 1) {
      cerr << label << ": " <<
	"WARNING! template named `" << (it->first)->getName() <<
	"' at address 0x" << it->first << " claims to have " <<
	it->second << "instances!" << endl;
    }

  }
}

/** List outstanding descriptors on stderr

    @param label (in) String to label output 
*/
void GA::GAServices::listDescriptors(std::string label) {
  std::map<DistArrayDescriptor *, int>::iterator it = descriptorRecord.begin();

  for ( ; it != descriptorRecord.end(); ++it) {
    if ( it->second == 1) {
      cerr << label << ": " <<
	"WARNING! descriptor `" << (it->first)->getName() <<
	"' at address 0x" << it->first << " still exists." << endl;
    }

    /** Report any bizarrities */

    else if ( it->second < 0 || it->second > 1) {
      cerr << label << ": " <<
	"WARNING! descriptor named `" << (it->first)->getName() <<
	"' at address 0x" << it->first << " claims to have " <<
	it->second << "instances!" << endl;
    }

  }
}

/** List outstanding arrays on stderr

    @param label (in) String to label output 
*/
void GA::GAServices::listArrays(std::string label) {
  std::map<DistArray *, int>::iterator it = arrayRecord.begin();

  for ( ; it != arrayRecord.end(); ++it) {
    if ( it->second == 1) {
      cerr << label << ": " <<
	"WARNING! array `" << (it->first)->getName() <<
	"' at address 0x" << it->first << " still exists." << endl;
    }

    /** Report any bizarrities */

    else if ( it->second < 0 || it->second > 1) {
      cerr << label << ": " <<
	"WARNING! array named `" << (it->first)->getName() <<
	"' at address 0x" << it->first << " claims to have " <<
	it->second << "instances!" << endl;
    }

  }
}

