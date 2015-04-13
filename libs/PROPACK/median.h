/*  Grassmann Averages
    Copyright (C) 2014  Søren Hauberg

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <cstdio>
#include <cstring>
#include <algorithm>
#include <vector>
#include <iostream>
//#include <tuple>

template<class value_type, class index_type>
value_type median (const value_type *data, index_type start, index_type stop)
{
  // Copy the data
//  std::cerr << __LINE__ << std::endl;
  std::vector<value_type> work(data+start, data+stop);
//  std::cerr << __LINE__ << std::endl;
  
  // Find median
  const index_type K = (stop - start) / 2;
  std::nth_element (work.begin(), work.begin()+K, work.end());
  const value_type median = work[K];
  
  return median;
}

template<class T1, class T2>
bool compare_first(std::pair<T1, T2> a, std::pair<T1, T2> b)
{
  return a.first < b.first;
}

template<class value_type, class index_type>
value_type weighted_median (const value_type *data, index_type start, index_type stop, const value_type *weights)
{
  // Copy the data
  typedef std::pair<value_type, value_type> weighted_value;
  const index_type N = stop - start;
  std::vector<weighted_value> work(N);
  value_type wsum = 0.0;
  for (index_type n = 0; n < N; n++)
    {
      work[n] = std::make_pair(data[start+n], weights[n]);
      wsum += weights[n];
    }
  wsum /= 2.0;
  
  // Find weighted median
  std::sort (work.begin(), work.end(), compare_first<value_type, value_type>);
  value_type accum = 0.0;
  index_type idx = 0;
  while (accum < wsum)
    accum += work[idx++].second;
  
  const value_type median = work[idx].first;
  
  return median;
}


