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

/*
#include <queue>

template<class value_type, class index_type>
value_type trimmed_mean (const value_type *data, index_type start, index_type stop, index_type K)
{
  // Some helpful typedefs
  typedef std::priority_queue< value_type, std::vector<value_type>, std::less<value_type> >    min_queue;
  typedef std::priority_queue< value_type, std::vector<value_type>, std::greater<value_type> > max_queue;
  
  // Step 0: Basic sanity checking; if we have less than 2K data elements return 0.
  const index_type numel = stop - start;
  if (numel < 2*K)
    return 0;
  
  // Step 1: find Kth smallest and largest elements of data
  max_queue max_q (data+start, data+start+K);
  min_queue min_q (data+start, data+start+K);
    
  for (index_type i = start+K; i < stop; i++)
    {
      const value_type data_i = data[i];
     
      min_q.push (data_i);
      min_q.pop ();
      
      max_q.push (data_i);
      max_q.pop ();
    }
  
  const value_type max_value = max_q.top ();
  const value_type min_value = min_q.top ();
  
  // Step 2: compute the mean of the inliers
  value_type N = 0;
  value_type mean = 0;
  for (index_type i = start; i < stop; i++)
    {
      const value_type data_i = data[i];
      if (data_i > min_value && data_i < max_value)
        {
          mean += data_i;
          N++;
        }
    }
  mean /= N;
  
  return mean;
}
*/

template<class value_type, class index_type>
value_type trimmed_mean (const value_type *data, index_type start, index_type stop, index_type K)
{
  // Step 0: Basic sanity checking; if we have less than 2K data elements return 0.
  const index_type numel = stop - start;
  if (numel < 2*K)
    return 0;

  // Copy the data
  std::vector<value_type> work(data+start, data+stop);
  
  // Step 1: find Kth smallest and largest elements of data
  std::nth_element (work.begin(), work.begin()+K, work.end(), std::less<value_type>());
  const value_type min_value = work[K];
  std::nth_element (work.begin(), work.begin()+K, work.end(), std::greater<value_type>());
  const value_type max_value = work[K];
  
  //std::cout << max_value << " " << min_value << std::endl;
  
  // Step 2: compute the mean of the inliers
  value_type N = 0;
  value_type mean = 0;
  for (index_type i = start; i < stop; i++)
    {
      const value_type data_i = data[i];
      if (data_i >= min_value && data_i <= max_value)
        {
          mean += data_i;
          N++;
        }
    }
  mean /= N;
  
  return mean;
}

