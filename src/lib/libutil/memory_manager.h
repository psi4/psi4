#ifndef _psi_src_bin_psimrccmemory_managerh_
#define _psi_src_bin_psimrccmemory_managerh_

#include <map>
#include <vector>
#include <string>


namespace psi{

/*
 * Computes the size in mebibytes (MiB) of a given amount of type T
 */
template <typename T>
double type_to_MiB(size_t n)
{
  // 1 MiB = 1048576 bytes
  size_t bites = n * static_cast<size_t>(sizeof(T));
  return(static_cast<double>(bites)/1048576.0);
}

/*
 * Convert bytes to mebibytes (MiB)
 */
double bytes_to_MiB(size_t n);

typedef struct {
	void*               variable;
	std::string         type;
	std::string         variableName;
	std::string         fileName;
	size_t              lineNumber;
	std::vector<size_t> argumentList;
} AllocationEntry;

class MemoryManager
{
public:
  MemoryManager(size_t maxcor = 256000000);
  ~MemoryManager();

  void MemCheck(FILE *output);

  size_t      get_FreeMemory()                   const {return(MaximumAllowed - CurrentAllocated);}
  size_t      get_CurrentAllocated()             const {return(CurrentAllocated);}
  size_t      get_MaximumAllowedMemory()         const {return(MaximumAllowed);}

  template <typename T>
  void allocate(const char *type, T*& matrix, size_t size, const char *variableName, const char *fileName, size_t lineNumber);
  template <typename T>
  void release_one(T*& matrix, const char *fileName, size_t lineNumber);

  template <typename T>
  void allocate(const char *type, T**& matrix, size_t size1, size_t size2, const char *variableName, const char *fileName, size_t lineNumber);
  template <typename T>
  void release_two(T**& matrix, const char *fileName, size_t lineNumber);

  template <typename T>
  void allocate(const char *type, T***& matrix,size_t size1,size_t size2,size_t size3, const char *variableName, const char *fileName, size_t lineNumber);
  template <typename T>
  void release_three(T***& matrix, const char *fileName, size_t lineNumber);
private:
  void RegisterMemory(void *mem, AllocationEntry& entry, size_t size);
  void UnregisterMemory(void *mem, size_t size, const char *fileName, size_t lineNumber);

  size_t CurrentAllocated;
  size_t MaximumAllocated;
  size_t MaximumAllowed;
  std::map<void *, AllocationEntry> AllocationTable;
};

template <typename T>
void MemoryManager::allocate(const char *type, T*& matrix, size_t size, const char *variableName, const char *fileName, size_t lineNumber)
{
  AllocationEntry newEntry;

  if(size<=0){
    matrix = NULL;
  }else{
    matrix    = new T[size];
    for(size_t i=0;i<size;i++)
      matrix[i]=static_cast<T>(0);   // Zero all the elements

    newEntry.variable = matrix;
    newEntry.type = type;
    newEntry.variableName = variableName;
    newEntry.fileName = fileName;
    newEntry.lineNumber = lineNumber;
    newEntry.argumentList.push_back(size);
    RegisterMemory(static_cast<void*>(matrix), newEntry, size*sizeof(T));
  }
}

template <typename T>
void MemoryManager::release_one(T*& matrix, const char *fileName, size_t lineNumber)
{
  if(matrix == NULL)
    return;

  size_t size = AllocationTable[static_cast<void*>(matrix)].argumentList[0];

  UnregisterMemory(static_cast<void*>(matrix), size*sizeof(T),fileName,lineNumber);

  delete[] matrix;
  matrix = NULL;
}

template <typename T>
void MemoryManager::allocate(const char *type, T**& matrix, size_t size1, size_t size2, const char *variableName, const char *fileName, size_t lineNumber)
{
  AllocationEntry newEntry;
  size_t size = size1*size2;

  if(size<=0){
    matrix = NULL;
    return;
  }else{
    matrix    = new T*[size1];
    T* vector = new T[size];
    for(size_t i=0;i<size;i++)
      vector[i]=static_cast<T>(0);   // Zero all the elements
    for(size_t i=0;i<size1;i++)
      matrix[i]=&(vector[i*size2]);  // Assign the rows pointers

    newEntry.variable = matrix;
    newEntry.type = type;
    newEntry.variableName = variableName;
    newEntry.fileName = fileName;
    newEntry.lineNumber = lineNumber;
    newEntry.argumentList.push_back(size1);
    newEntry.argumentList.push_back(size2);
    RegisterMemory(static_cast<void*>(matrix), newEntry, size*sizeof(T));
  }
}

template <typename T>
void MemoryManager::release_two(T**& matrix, const char *fileName, size_t lineNumber)
{
  if(matrix == NULL)
    return;

  size_t size = AllocationTable[static_cast<void*>(matrix)].argumentList[0] * AllocationTable[static_cast<void*>(matrix)].argumentList[1];

  UnregisterMemory(static_cast<void*>(matrix), size*sizeof(T),fileName,lineNumber);

  delete[] matrix[0];
  delete[] matrix;
  matrix = NULL;
}

template <typename T>
void MemoryManager::allocate(const char *type, T***& matrix,size_t size1,size_t size2,size_t size3, const char *variableName, const char *fileName, size_t lineNumber)
{
  AllocationEntry newEntry;
  size_t size = size1*size2*size3;
  if(size<=0){
    matrix = NULL;
    return;
  }else{
    matrix    = new T**[size1];
    for(size_t i=0;i<size1;i++)
      matrix[i]= new T*[size2];
    T* vector = new T[size];
    for(size_t i=0;i<size;i++)
      vector[i]=static_cast<T>(0);   // Zero all the elements
    for(size_t i=0;i<size1;i++)
      for(size_t j=0;j<size2;j++)
        matrix[i][j]=&(vector[i*size2*size3+j*size3]);  // Assign the rows pointers
    newEntry.variable = matrix;
    newEntry.type = type;
    newEntry.variableName = variableName;
    newEntry.fileName = fileName;
    newEntry.lineNumber = lineNumber;
    newEntry.argumentList.push_back(size1);
    newEntry.argumentList.push_back(size2);
    newEntry.argumentList.push_back(size3);
    RegisterMemory(static_cast<void*>(matrix), newEntry, size*sizeof(T));
  }
}

template <typename T>
void MemoryManager::release_three(T***& matrix, const char *fileName, size_t lineNumber)
{
  if(matrix == NULL)
    return;

  size_t size1 = AllocationTable[static_cast<void*>(matrix)].argumentList[0];
  size_t size = AllocationTable[static_cast<void*>(matrix)].argumentList[0] * AllocationTable[static_cast<void*>(matrix)].argumentList[1]
              * AllocationTable[static_cast<void*>(matrix)].argumentList[2];

  UnregisterMemory(static_cast<void*>(matrix), size*sizeof(T),fileName,lineNumber);

  delete[] matrix[0][0];
  for(size_t i=0;i<size1;i++)
    delete[] matrix[i];
  delete[] matrix;
  matrix = NULL;
}

#define allocate1(type, variable, size) \
  memory_manager->allocate(#type, variable, size, #variable, __FILE__, __LINE__);
#define release1(variable) \
  memory_manager->release_one(variable, __FILE__, __LINE__);

#define allocate2(type, variable, size1, size2) \
  memory_manager->allocate(#type, variable, size1, size2, #variable, __FILE__, __LINE__);
#define release2(variable) \
  memory_manager->release_two(variable, __FILE__, __LINE__);

#define allocate3(type, variable, size1, size2, size3) \
  memory_manager->allocate(#type, variable, size1, size2, size3, #variable, __FILE__, __LINE__);
#define release3(variable) \
  memory_manager->release_three(variable, __FILE__, __LINE__);

} /* End Namespaces */

#endif // _psi_src_bin_psimrccmemory_managerh_
