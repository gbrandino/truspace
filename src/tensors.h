/** @file tensors.h 
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
 * \date October 2017
 *          
 */
#ifndef TENSORS_H
#define TENSORS_H
#include <iostream>
#include <cstdlib>

/** Template class for 2D tensors
 * This classes do not provide an allocator. They are allocated automatically using the costum constructor tesor?D( int int ....)
 * The default constructor can be used, but only to create dummy tensors, that needs to be assigned using the assignment operator (=)
 * 
 * Memebers are all public.  
 *   - n1, n2, ... contains the dimensions
 *   - the pointer data contains the data in col major order, such that it can be pass to BLAS/LAPACK easily as A.data
 *   - reading and writing single elements is done using the () operator so element 1,2,5 of tensor A is accessed A(1,2,5) 
 */
template <class T>
class tensor2D
{
    public:
        tensor2D();
        tensor2D(int, int);
        ~tensor2D();
        tensor2D<T>& operator=(const tensor2D<T>&);
        T& operator()(const int,const int);
        T operator()(const int,const int) const;
        /** First ensor size*/
        int n1; 
        /** Second size*/
        int n2; 
        /** Data segment */
        T* data; 
        
};

/** Default Constructor, it allocates a 1x1 tensor 

 */
template <class T> tensor2D<T>::tensor2D()
{
    data = new T[1];
}

/** Constructor to allocate custom size 2D tensor.
 * @param[in] n1_in First dimension
 * @param[in] n2_in Second dimension
 */
template <class T>
tensor2D<T>::tensor2D(int n1_in, int n2_in)
:
n1(n1_in),
n2(n2_in)
{
    data = new T[n1*n2];
}

/** Destructor, deallocates data member.
 */
template <class T>
tensor2D<T>::~tensor2D()
{
    delete [] data;
}

/** Assignment operator.
 * @param[in] tens The tensor to copy
 * return A reference to the current object, to allow expressions such as a=b=c to work
 */
template<class T>
tensor2D<T>& tensor2D<T>::operator=(const tensor2D<T>& tens)
{
    if (this != &tens)
        {
        delete [] data;
        n1=tens.n1;
        n2=tens.n2;
        int size = n1*n2;
        data=new T[size];
            for (int i=0; i<size;i++)
                            data[i]=tens.data[i];
        }
    return *this;   
}

/** Element access operator, with bound check.
 * @param[in] i1 First index of the element to access
 * @param[in] i2 Second index of the element to access
 */ 
template <class T>
T& tensor2D<T>::operator()(const int i1, const int i2)
{
    if (i1 > (n1-1) or i2 > (n2-1) )
        {
        std::cerr << " ERROR: Out of bound on 2D tensor! Requested " << std::endl;
        std::cerr << i1 << " (size " << n1 << ")" << std::endl;
        std::cerr << i2 << " (size " << n2 << ")" << std::endl;
        std::exit(1);
        }
    else
        {
        return data[i1*n2 + i2];
        } 
}

/** Element access operator, with bound check - Constant tensor version.
 * @param[in] i1 First index of the element to access
 * @param[in] i2 Second index of the element to access
 */
template <class T>
T tensor2D<T>::operator()(const int i1, const int i2) const
{
    if (i1 > (n1-1) or i2 > (n2-1) )
        {
        std::cerr << " ERROR: Out of bound on 2D tensor! Requested " << std::endl;
        std::cerr << i1 << " (size " << n1 << ")" << std::endl;
        std::cerr << i2 << " (size " << n2 << ")" << std::endl;
        std::exit(1);
        }
    else
        {
        return data[i1*n2 + i2];
        }
}

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
class tensor3D
{
    public:
        tensor3D();
        tensor3D(int, int, int);
        ~tensor3D();
        tensor3D<T>& operator=(const tensor3D<T>&);
        T& operator()(const int,const int, const int);
        T operator()(const int,const int, const int) const;
        /** First tensor size*/
        int n1; 
        /** Second tensor size*/
        int n2; 
        /** Third tensor size*/
        int n3; 
        /** Data segment */
        T* data; 

};

/** Default Constructor, it allocates a 1x1x1 tensor 
 */
template <class T>
tensor3D<T>::tensor3D()
{
    data = new T[1];
}

/** Constructor to allocate custom size 3D tensor.
 * @param[in] n1_in First dimension
 * @param[in] n2_in Second dimension
 * @param[in] n3_in Thrid dimension
 */
template <class T>
tensor3D<T>::tensor3D(int n1_in, int n2_in, int n3_in)
:
n1(n1_in),
n2(n2_in),
n3(n3_in)
{
    data = new T[n1*n2*n3];
}

/** Destructor, deallocates data member.
 */
template <class T>
tensor3D<T>::~tensor3D()
{
    delete [] data;
}

/** Assignment operator.
 * @param[in] tens The tensor to copy
 * return A reference to the current object, to allow expressions such as a=b=c to work
 */
template<class T>
tensor3D<T>& tensor3D<T>::operator=(const tensor3D<T>& tens)
{
    if (this != &tens)
        {
        delete [] data;
        n1=tens.n1;
        n2=tens.n2;
        n3=tens.n3;
        int size = n1*n2*n3;
        data=new T[size];
            for (int i=0; i<size;i++)
                            data[i]=tens.data[i];
        }
    return *this;
}


/** Element access operator, with bound check.
 * @param[in] i1 First index of the element to access
 * @param[in] i2 Second index of the element to access
 * @param[in] i3 Third index of the element to access
 */
template <class T>
T& tensor3D<T>::operator()(const int i1, const int i2, const int i3)
{
    if (i1 > (n1-1) or i2 > (n2-1) or i3 > (n3-1))
        {
        std::cerr << " ERROR: Out of bound on 3D tensor! Requested " << std::endl;
        std::cerr << i1 << " (size " << n1 << ")" << std::endl;
        std::cerr << i2 << " (size " << n2 << ")" << std::endl;
        std::cerr << i3 << " (size " << n3 << ")" << std::endl;
        std::exit(1);
        }
    else
        {
        return data[i1*n2*n3 + i2*n3 + i3 ];
        }
}

/** Element access operator, with bound check - Constant tensor version.
 * @param[in] i1 First index of the element to access
 * @param[in] i2 Second index of the element to access
 * @param[in] i3 Third index of the element to access
 */
template <class T>
T tensor3D<T>::operator()(const int i1, const int i2, const int i3) const
{
    if (i1 > (n1-1) or i2 > (n2-1) or i3 > (n3-1))
        {
        std::cerr << " ERROR: Out of bound on 3D tensor! Requested " << std::endl;
        std::cerr << i1 << " (size " << n1 << ")" << std::endl;
        std::cerr << i2 << " (size " << n2 << ")" << std::endl;
        std::cerr << i3 << " (size " << n3 << ")" << std::endl;
        std::exit(1);
        }
    else
        {
        return data[i1*n2*n3 + i2*n3 + i3 ];
        }
}


/** Template class for 4D tensors
 * This classes do not provide an allocator. They are allocated automatically using the costum constructor tesor?D( int int ....)
 * The default constructor can be used, but only to create dummy tensors, that needs to be assigned using the assignment operator (=)
 * 
 * Memebers are all public.  
 *   - n1, n2, ... contains the dimensions
 *   - the pointer data contains the data in col major order, such that it can be pass to BLAS/LAPACK easily as A.data
 *   - reading and writing single elements is done using the () operator so element 1,2,5 of tensor A is accessed A(1,2,5)
 */
template <class T>
class tensor4D
{
    public:
        tensor4D();
        tensor4D(int, int, int, int);
        ~tensor4D();
        tensor4D<T>& operator=(const tensor4D<T>&);
        T& operator()(const int,const int, const int, const int);
        T operator()(const int,const int, const int, const int) const;
        /** First tensor size*/
        int n1; 
        /** Second tensor size*/
        int n2; 
        /** Third tensor size*/
        int n3; 
        /** Fourth tensor size */
        int n4;
        /** Data segment */
        T* data;

};

/** Default Constructor, it allocates a 1x1x1x1 tensor 
 */
template <class T>
tensor4D<T>::tensor4D()
{
    data = new T[1];
}


/** Constructor to allocate custom size 4D tensor.
 * @param[in] n1_in First dimension
 * @param[in] n2_in Second dimension
 * @param[in] n3_in Third dimension
 * @param[in] n4_in Fourth dimension
 */
template <class T>
tensor4D<T>::tensor4D(int n1_in, int n2_in, int n3_in, int n4_in)
:
n1(n1_in),
n2(n2_in),
n3(n3_in),
n4(n4_in)
{
    data = new T[n1*n2*n3*n4];
}

/** Destructor, deallocates data member.
 */
template <class T>
tensor4D<T>::~tensor4D()
{
    delete [] data;
}

/** Assignment operator.
 * @param[in] tens The tensor to copy
 * return A reference to the current object, to allow expressions such as a=b=c to work
 */
template<class T>
tensor4D<T>& tensor4D<T>::operator=(const tensor4D<T>& tens)
{
    if (this != &tens)
        {
        delete [] data;
        n1=tens.n1;
        n2=tens.n2;
        n3=tens.n3;
        n4=tens.n4;
        int size = n1*n2*n3*n4;
        data=new T[size];
            for (int i=0; i<size;i++)
                            data[i]=tens.data[i];
        }
    return *this;
}


/** Element access operator, with bound check.
 * @param[in] i1 First index of the element to access
 * @param[in] i2 Second index of the element to access
 * @param[in] i3 Third index of the element to access
 * @param[in] i4 Fourth index of the element to access
 */
template <class T>
T& tensor4D<T>::operator()(const int i1, const int i2, const int i3, const int i4)
{
    if (i1 > (n1-1) or i2 > (n2-1) or i3 > (n3-1) or i4 > (n4-1))
        {
        std::cerr << " ERROR: Out of bound on 4D tensor! Requested " << std::endl;
        std::cerr << i1 << " (size " << n1 << ")" << std::endl;
        std::cerr << i2 << " (size " << n2 << ")" << std::endl;
        std::cerr << i3 << " (size " << n3 << ")" << std::endl;
        std::cerr << i4 << " (size " << n4 << ")" << std::endl;
        std::exit(1);
        }
    else
        {
        return data[i1*n2*n3*n4 + i2*n3*n4 + i3*n4 + i4 ];
        }
}

/** Element access operator, with bound check - Constant tensor version.
 * @param[in] i1 First index of the element to access
 * @param[in] i2 Second index of the element to access
 * @param[in] i3 Third index of the element to access
 * @param[in] i4 Fourth index of the element to access
 */
template <class T>
T tensor4D<T>::operator()(const int i1, const int i2, const int i3, const int i4) const
{
    if (i1 > (n1-1) or i2 > (n2-1) or i3 > (n3-1) or i4 > (n4-1))
        {
        std::cerr << " ERROR: Out of bound on 4D tensor! Requested " << std::endl;
        std::cerr << i1 << " (size " << n1 << ")" << std::endl;
        std::cerr << i2 << " (size " << n2 << ")" << std::endl;
        std::cerr << i3 << " (size " << n3 << ")" << std::endl;
        std::cerr << i4 << " (size " << n4 << ")" << std::endl;
        std::exit(1);
        }
    else
        {
        return data[i1*n2*n3*n4 + i2*n3*n4 + i3*n4 + i4 ];
        }
}

/** Template class for 5D tensors
 * This classes do not provide an allocator. They are allocated automatically using the costum constructor tesor?D( int int ....)
 * The default constructor can be used, but only to create dummy tensors, that needs to be assigned using the assignment operator (=)
 * 
 * Memebers are all public.  
 *   - n1, n2, ... contains the dimensions
 *   - the pointer data contains the data in col major order, such that it can be pass to BLAS/LAPACK easily as A.data
 *   - reading and writing single elements is done using the () operator so element 1,2,5 of tensor A is accessed A(1,2,5) 
 */
template <class T>
class tensor5D
{
    public:
        tensor5D();
        tensor5D(int, int, int, int, int);
        ~tensor5D(); 
        tensor5D<T>& operator=(const tensor5D<T>&);
        T& operator()(const int,const int, const int, const int, const int);
        T operator()(const int,const int, const int, const int, const int) const;
        /** First tensor size*/
        int n1;
        /** Second tensor size*/
        int n2;
        /** Third tensor size*/
        int n3;
        /** Fourth tensor size */
        int n4;
        /** Fifth tensor size */
        int n5;
        /** Data segment */
        T* data;
};

/** Default Constructor, it allocates a 1x1x1x1x1 tensor 
 */
template <class T>
tensor5D<T>::tensor5D()
{
    data = new T[1];
}

/** Constructor to allocate custom size 5D tensor.
 * @param[in] n1_in First dimension
 * @param[in] n2_in Second dimension
 * @param[in] n3_in Third dimension
 * @param[in] n4_in Fourth dimension
 * @param[in] n5_in Fifth dimension
 */
template <class T>
tensor5D<T>::tensor5D(int n1_in, int n2_in, int n3_in, int n4_in, int n5_in)
:
n1(n1_in),
n2(n2_in),
n3(n3_in),
n4(n4_in),
n5(n5_in)
{   
    data = new T[n1*n2*n3*n4*n5];
}

/** Destructor, deallocates data member.
 */
template <class T>
tensor5D<T>::~tensor5D()
{
    delete [] data;
}


/** Assignment operator.
 * @param[in] tens The tensor to copy
 * return A reference to the current object, to allow expressions such as a=b=c to work
 */
template<class T>
tensor5D<T>& tensor5D<T>::operator=(const tensor5D<T>& tens)
{
    if (this != &tens)
        {
        delete [] data;
        n1=tens.n1;
        n2=tens.n2;
        n3=tens.n3;
        n4=tens.n4;
        n5=tens.n5;
        int size = n1*n2*n3*n4*n5;
        data=new T[size];
            for (int i=0; i<size;i++)
                            data[i]=tens.data[i];
        }
    return *this;
}

/** Element access operator, with bound check. 
 * @param[in] i1 First index of the element to access
 * @param[in] i2 Second index of the element to access
 * @param[in] i3 Third index of the element to access
 * @param[in] i4 Fourth index of the element to access
 * @param[in] i5 Fifth index of the element to access
 */
template <class T>
T& tensor5D<T>::operator()(const int i1, const int i2, const int i3, const int i4, const int i5)
{
    if (i1 > (n1-1) or i2 > (n2-1) or i3 > (n3-1) or i4 > (n4-1) or i5 > (n5-1) )
        {
        std::cerr << " ERROR: Out of bound on 5D tensor! Requested " << std::endl;
        std::cerr << i1 << " (size " << n1 << ")" << std::endl;
        std::cerr << i2 << " (size " << n2 << ")" << std::endl;
        std::cerr << i3 << " (size " << n3 << ")" << std::endl;
        std::cerr << i4 << " (size " << n4 << ")" << std::endl;
        std::cerr << i5 << " (size " << n5 << ")" << std::endl;
        std::exit(1);
        }
    else
        {
        return data[i1*n2*n3*n4*n5 + i2*n3*n4*n5 + i3*n4*n5 + i4*n5 + i5];
        }
}

/** Element access operator, with bound check - Constant tensor version.
 * @param[in] i1 First index of the element to access
 * @param[in] i2 Second index of the element to access 
 * @param[in] i3 Third index of the element to access
 * @param[in] i4 Fourth index of the element to access
 * @param[in] i5 Fifth index of the element to access
 */
template <class T>
T tensor5D<T>::operator()(const int i1, const int i2, const int i3, const int i4, const int i5) const
{
    if (i1 > (n1-1) or i2 > (n2-1) or i3 > (n3-1) or i4 > (n4-1) or i5 > (n5-1) )
        {
        std::cerr << " ERROR: Out of bound on 5D tensor! Requested " << std::endl;
        std::cerr << i1 << " (size " << n1 << ")" << std::endl;
        std::cerr << i2 << " (size " << n2 << ")" << std::endl;
        std::cerr << i3 << " (size " << n3 << ")" << std::endl;
        std::cerr << i4 << " (size " << n4 << ")" << std::endl;
        std::cerr << i5 << " (size " << n5 << ")" << std::endl;
        std::exit(1);
        }
    else
        {
        return data[i1*n2*n3*n4*n5 + i2*n3*n4*n5 + i3*n4*n5 + i4*n5 + i5];
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
class tensor6D
{
    public:
        tensor6D();
        tensor6D(int, int, int, int, int, int);
        ~tensor6D();
        tensor6D<T>& operator=(const tensor6D<T>&);
        T& operator()(const int,const int, const int, const int, const int, const int);
        T operator()(const int,const int, const int, const int, const int, const int) const;
        /** First tensor size*/
        int n1;
        /** Second tensor size*/
        int n2;
        /** Third tensor size*/
        int n3;
        /** Fourth tensor size */
        int n4;
        /** Fifth tensor size */
        int n5;
        /** Sixth tensor size */
        int n6; 
        /** Data segment */ 
        T* data;

};

/** Default Constructor, it allocates a 1x1x1x1x1x1 tensor 
 */
template <class T>
tensor6D<T>::tensor6D()
{
    data = new T[1];
}

/** Constructor to allocate custom size 2D tensor.
 * @param[in] n1_in First dimension
 * @param[in] n2_in Second dimension
 * @param[in] n3_in Third dimension
 * @param[in] n4_in Fourth dimension
 * @param[in] n5_in Fifth dimension
 * @param[in] n6_in Sixth dimension
 */
template <class T>
tensor6D<T>::tensor6D(int n1_in, int n2_in, int n3_in, int n4_in, int n5_in, int n6_in)
:
n1(n1_in),
n2(n2_in),
n3(n3_in),
n4(n4_in),
n5(n5_in),
n6(n6_in)
{
    data = new T[n1*n2*n3*n4*n5*n6];
}

/** Destructor, deallocates data member.
 */
template <class T>
tensor6D<T>::~tensor6D()
{
    delete [] data;
}

/** Assignment operator.
 * @param[in] tens The tensor to copy
 * return A reference to the current object, to allow expressions such as a=b=c to work
 */
template<class T>
tensor6D<T>& tensor6D<T>::operator=(const tensor6D<T>& tens)
{
    if (this != &tens)
        {
        delete [] data;
        n1=tens.n1;
        n2=tens.n2;
        n3=tens.n3;
        n4=tens.n4;
        n5=tens.n5;
        n6=tens.n6;
        int size = n1*n2*n3*n4*n5*n6;
        data=new T[size];
            for (int i=0; i<size;i++)
                            data[i]=tens.data[i];
        }
    return *this;
}


/** Element access operator, with bound check. 
 * @param[in] i1 First index of the element to access
 * @param[in] i2 Second index of the element to access 
 * @param[in] i3 Third index of the element to access
 * @param[in] i4 Fourth index of the element to access
 * @param[in] i5 Fifth index of the element to access
 * @param[in] i6 Sixth index of the element to access
 */
template <class T>
T& tensor6D<T>::operator()(const int i1, const int i2, const int i3, const int i4, const int i5, const int i6)
{
    if (i1 > (n1-1) or i2 > (n2-1) or i3 > (n3-1) or i4 > (n4-1) or i5 > (n5-1) or i6 > (n6-1))
        {
        std::cerr << " ERROR: Out of bound on 6D tensor! Requested " << std::endl;
        std::cerr << i1 << " (size " << n1 << ")" << std::endl;
        std::cerr << i2 << " (size " << n2 << ")" << std::endl;
        std::cerr << i3 << " (size " << n3 << ")" << std::endl;
        std::cerr << i4 << " (size " << n4 << ")" << std::endl;
        std::cerr << i5 << " (size " << n5 << ")" << std::endl;
        std::cerr << i6 << " (size " << n6 << ")" << std::endl;
        std::exit(1);
        }
    else
        {
        return data[i1*n2*n3*n4*n5*n6 + i2*n3*n4*n5*n6 + i3*n4*n5*n6 + i4*n5*n6 + i5*n6 + i6 ];
        }
}

/** Element access operator, with bound check - Constant tensor version.
 * @param[in] i1 First index of the element to access
 * @param[in] i2 Second index of the element to access 
 * @param[in] i3 Third index of the element to access
 * @param[in] i4 Fourth index of the element to access
 * @param[in] i5 Fifth index of the element to access
 * @param[in] i6 Sixth index of the element to access
 */
template <class T>
T tensor6D<T>::operator()(const int i1, const int i2, const int i3, const int i4, const int i5, const int i6) const
{
    if (i1 > (n1-1) or i2 > (n2-1) or i3 > (n3-1) or i4 > (n4-1) or i5 > (n5-1) or i6 > (n6-1))
        {
        std::cerr << " ERROR: Out of bound on 6D tensor! Requested " << std::endl;
        std::cerr << i1 << " (size " << n1 << ")" << std::endl;
        std::cerr << i2 << " (size " << n2 << ")" << std::endl;
        std::cerr << i3 << " (size " << n3 << ")" << std::endl;
        std::cerr << i4 << " (size " << n4 << ")" << std::endl;
        std::cerr << i5 << " (size " << n5 << ")" << std::endl;
        std::cerr << i6 << " (size " << n6 << ")" << std::endl;
        std::exit(1);
        }
    else
        {
        return data[i1*n2*n3*n4*n5*n6 + i2*n3*n4*n5*n6 + i3*n4*n5*n6 + i4*n5*n6 + i5*n6 + i6 ];
        }
}

#endif


