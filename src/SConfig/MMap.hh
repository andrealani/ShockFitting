// Copyright (C) 2010 Andrea Lani
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef MMap_hh
#define MMap_hh

#include <vector>
#include <cassert>
#include <numeric>
#include <algorithm>

//----------------------------------------------------------------------------//

namespace SConfig {
  
//----------------------------------------------------------------------------//

// This class implements a more efficient std::map based on a std::vector
// provided with sorting and search operations
// @author Andrea Lani
template <typename KEY, typename VALUE>
class MMap {
public:
  
  // constructor
  MMap(size_t NKEYS = 0) : m_usable(false), m_map()
  {
    reserve(NKEYS);
  }
  
  // destructor
  ~MMap() {}
  
  // Insert pair <KEY, VALUE>
  void insert(const KEY& key, const VALUE& value) 
  {    
    m_map.push_back(std::pair<KEY, VALUE>(key, value)); m_usable = false;
  }
  
  // find the given key: if available keyFound is set to "true"
  typedef typename std::vector<std::pair<KEY, VALUE> >::iterator MapItr;
  std::pair<MapItr, MapItr> find(const KEY& key, bool& keyFound)
  {
    if (!m_usable) {sort();} // first sort if needed
    
    keyFound = false;
    std::pair<MapItr, MapItr> mapValue;
    if (m_map.empty()) {return mapValue;}
    
    if (std::binary_search (m_map.begin(),m_map.end(), key, Compare())) {
      keyFound = true;
      mapValue = std::equal_range(m_map.begin(),m_map.end(), key, Compare());
      // keyFound = (mapValue.first->first == key);
    }
    return mapValue;
  }
  
  // get the total numbers of keys
  size_t size() const {return m_map.size();}
  
  // preallocate the memory for the map
  void reserve(size_t nkeys)
  {
    if (nkeys > 0) m_map.reserve(nkeys);
  }
  
  // erase all entries in the map
  void clear()
  {
    std::vector<std::pair<KEY,VALUE> >().swap(m_map);
  }
  
  // sort all the keys (this is needed before using find())
  void sort() 
  { 
    std::sort(m_map.begin(),m_map.end(), LessThan());
    m_usable = true;
  }
  
  // get the pair corresponding to the given position
  // with no guarantee that the map has been sorted
  std::pair<KEY,VALUE> getPair(size_t pos) const {return m_map[pos];}
  
  /// overloading of the operator"[]" for assignment
  VALUE& operator[] (size_t pos) {
    assert(pos < size()); 
    return m_map[pos].second;
  }
  
  // map iterator pointing at the first entry
  MapItr begin() {return m_map.begin();}
  
  // map iterator pointing at the last entry
  MapItr end() {return m_map.end();}
  
private: 
  
  // This functor class allows to sort two given pairs so that their
  // keys are in increasing order.
  /// @author Andrea Lani
  class LessThan {
  public:
    //	overloading of the operator() for making this class a functor
    bool operator() (const std::pair<KEY,VALUE>& a,
		     const std::pair<KEY,VALUE>& b) const
    {
      return (a.first < b.first) ? true : false;
    }
  };
  
  // this functor class compares a pair and a key and viceversa
  class Compare {
  public:
    
    //	overloading of the operator() for making this class a functor
    bool operator() (const std::pair<KEY,VALUE>& p1, const KEY& key) const
    {
      return (p1.first < key) ? true : false;
    }
    
    //	overloading of the operator() for making this class a functor
    bool operator() (const KEY& key, const std::pair<KEY,VALUE>& p1) const
    {
      return (p1.first > key) ? true : false;
    }
  };
  
private: //data
  
  // check if the map is usable (meaning already sorted)
  bool m_usable;
  
  // array storing all pairs<KEY, VALUE>
  std::vector<std::pair<KEY, VALUE> > m_map;
  
};

//----------------------------------------------------------------------------//

} // end namespace SConfig

//----------------------------------------------------------------------------//

#endif
