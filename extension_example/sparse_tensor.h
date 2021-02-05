/** @file sparse_tensors.h 
 * \brief Set of helper templates class to manage 2D, 3D, 4D, 5D, 6D tensors 
 * 
 * This classes do not provide an allocator. They are allocated automatically using the costum constructor tesor?D( int int ....)
 * The default constructor can be used, but only to create dummy tensors, that needs to be assigned using the assignment operator (=)
 * 
 * Memebers are all public.  
 *   - n1, n2, ... contains the dimensions
 *   - the pointer data contains the data in col major order, such that it can be pass to BLAS/LAPACK easily as A.data
 *   - reading and writing single elements is done using the () operator so element 1,2,5 of tensor A is accessed A(1,2,5)
 * 
 * \author Giuseppe Piero Brandino - eXact-lab s.r.l.
 * \date May 2019
 *          
 */
#ifndef SPARSE_TENSORS_H
#define SPARSE_TENSORS_H
#include <iostream>
#include <cstdlib>
#include <unordered_map>
#include <tuple>
#include <limits>

typedef std::tuple<int, int, int> tuple3_key_t;

struct key_hash3 : public std::unary_function<tuple3_key_t, std::size_t>
{
 std::size_t operator()(const tuple3_key_t& k) const
 {
   return std::get<0>(k) ^ std::get<1>(k) ^ std::get<2>(k);
 }
};

typedef std::tuple<int, int, int, int, int, int> tuple6_key_t;

struct key_hash6 : public std::unary_function<tuple6_key_t, std::size_t>
{
 std::size_t operator()(const tuple6_key_t& k) const
 {
   return std::get<0>(k) ^ std::get<1>(k) ^ std::get<2>(k) ^ std::get<3>(k) ^ std::get<4>(k) ^ std::get<5>(k) ;
 }
};

/** Template class for 3D tensors
 * This classes do not provide an allocator. They are allocated automatically using the costum constructor tesor?D( int int ....)
 * The default constructor can be used, but only to create dummy tensors, that needs to be assigned using the assignment operator (=)
 * 
 * Memebers are all public.  
 *   - n1, n2, ... contains the dimensions
 *   - the pointer data contains the data in col major order, such that it can be pass to BLAS/LAPACK easily as A.data
 *   - reading and writing single elements is done using the () operator so element 1,2,5 of tensor A is accessed A(1,2,5) 
 */
template <class T>
class sparse_tensor3D
{
    public:
        sparse_tensor3D();
        ~sparse_tensor3D(); 
        sparse_tensor3D<T>& operator=(const sparse_tensor3D<T>&);
        T& set(const int,const int, const int);
        T operator()(const int,const int, const int);
        T operator()(const int,const int, const int) const;
        std::unordered_map< tuple3_key_t, T, key_hash3 > data;
        int size;

};

/** Default Constructor, it allocates a 1x1x1x1x1x1 tensor 
 */
template <class T>
sparse_tensor3D<T>::sparse_tensor3D()
{
    size = 0;
}


/** Destructor, deallocates data member.
 */
template <class T>
sparse_tensor3D<T>::~sparse_tensor3D()
{
}

/** Assignment operator.
 * @param[in] tens The tensor to copy
 * return A reference to the current object, to allow expressions such as a=b=c to work
 */
template<class T>
sparse_tensor3D<T>& sparse_tensor3D<T>::operator=(const sparse_tensor3D<T>& tens)
{
    if (this != &tens)
        {
            data = tens.data;
            size = tens.size;
        }
    return *this;
}

template <class T>
T& sparse_tensor3D<T>::set(const int i1, const int i2, const int i3)
{
    if (data.find(std::make_tuple(i1, i2, i3)) == data.end())
        size++;
    return data[std::make_tuple(i1, i2, i3)];
}




/** getter operator, with presence check. 
 * @param[in] i1 First index of the element to access
 * @param[in] i2 Second index of the element to access 
 * @param[in] i3 Third index of the element to access
 */
template <class T>
T sparse_tensor3D<T>::operator()(const int i1, const int i2, const int i3)
{
    if (data.find(std::make_tuple(i1, i2, i3)) != data.end())
        return data[std::make_tuple(i1, i2, i3)];
    else{
        std::cerr << "ERROR Key ( " << i1 << ", " << i2 << ", " << i3 << " ) not present! Returning NaN" << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }
}

/** getter operator, with bound check - Constant tensor version.
 * @param[in] i1 First index of the element to access
 * @param[in] i2 Second index of the element to access 
 * @param[in] i3 Third index of the element to access
 */
template <class T>
T sparse_tensor3D<T>::operator()(const int i1, const int i2, const int i3) const
{
    if (data.find(std::make_tuple(i1, i2, i3)) != data.end())
        return data[std::make_tuple(i1, i2, i3)];
    else{
        std::cerr << "ERROR: Key ( " << i1 << ", " << i2 << ", " << i3 << " ) not present! Returning to NaN" << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }
}



/** Template class for 6D tensors
 * This classes do not provide an allocator. They are allocated automatically using the costum constructor tesor?D( int int ....)
 * The default constructor can be used, but only to create dummy tensors, that needs to be assigned using the assignment operator (=)
 * 
 * Memebers are all public.  
 *   - n1, n2, ... contains the dimensions
 *   - the pointer data contains the data in col major order, such that it can be pass to BLAS/LAPACK easily as A.data
 *   - reading and writing single elements is done using the () operator so element 1,2,5 of tensor A is accessed A(1,2,5) 
 */
template <class T>
class sparse_tensor6D
{
    public:
        sparse_tensor6D();
        ~sparse_tensor6D();
        sparse_tensor6D<T>& operator=(const sparse_tensor6D<T>&);
        T& set(const int,const int, const int, const int, const int, const int);
        T operator()(const int,const int, const int, const int, const int, const int);
        T operator()(const int,const int, const int, const int, const int, const int) const;
        std::unordered_map< tuple6_key_t, T, key_hash6 > data;
        int size;

};

/** Default Constructor, it allocates a 1x1x1x1x1x1 tensor 
 */
template <class T>
sparse_tensor6D<T>::sparse_tensor6D()
{
    size = 0;
}


/** Destructor, deallocates data member.
 */
template <class T>
sparse_tensor6D<T>::~sparse_tensor6D()
{
}

/** Assignment operator.
 * @param[in] tens The tensor to copy
 * return A reference to the current object, to allow expressions such as a=b=c to work
 */
template<class T>
sparse_tensor6D<T>& sparse_tensor6D<T>::operator=(const sparse_tensor6D<T>& tens)
{
    if (this != &tens)
        {
            data = tens.data;
            size = tens.size;
        }
    return *this;
}

template <class T>
T& sparse_tensor6D<T>::set(const int i1, const int i2, const int i3, const int i4, const int i5, const int i6)
{
    if (data.find(std::make_tuple(i1, i2, i3, i4, i5, i6)) == data.end())
        size++;
    return data[std::make_tuple(i1, i2, i3, i4, i5, i6)];
}




/** getter operator, with presence check. 
 * @param[in] i1 First index of the element to access
 * @param[in] i2 Second index of the element to access 
 * @param[in] i3 Third index of the element to access
 * @param[in] i4 Fourth index of the element to access
 * @param[in] i5 Fifth index of the element to access
 * @param[in] i6 Sixth index of the element to access
 */
template <class T>
T sparse_tensor6D<T>::operator()(const int i1, const int i2, const int i3, const int i4, const int i5, const int i6)
{
    if (data.find(std::make_tuple(i1, i2, i3, i4, i5, i6)) != data.end()) 
        return data[std::make_tuple(i1, i2, i3, i4, i5, i6)];
    else{
        std::cerr << "ERROR Key ( " << i1 << ", " << i2 << ", " << i3 << ", " << i4 << ", " << i5 << ", " << i6 << " ) not present! Returning NaN" << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }
}

/** getter operator, with bound check - Constant tensor version.
 * @param[in] i1 First index of the element to access
 * @param[in] i2 Second index of the element to access 
 * @param[in] i3 Third index of the element to access
 * @param[in] i4 Fourth index of the element to access
 * @param[in] i5 Fifth index of the element to access
 * @param[in] i6 Sixth index of the element to access
 */
template <class T>
T sparse_tensor6D<T>::operator()(const int i1, const int i2, const int i3, const int i4, const int i5, const int i6) const
{
    if (data.find(std::make_tuple(i1, i2, i3, i4, i5, i6)) != data.end())    
        return data[std::make_tuple(i1, i2, i3, i4, i5, i6)];
    else{
        std::cerr << "ERROR: Key ( " << i1 << ", " << i2 << ", " << i3 << ", " << i4 << ", " << i5 << ", " << i6 << " ) not present! Returning to NaN" << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }
}

#endif


