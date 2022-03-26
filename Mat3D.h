// W. H. Press, et al, "Numerical Recipes"

#ifndef __Mat3D_h__
#define __Mat3D_h__

template <class T>
class Mat3d {
private:
    int nn;
    int mm;
    int kk;
    T ***v;
public:
    Mat3d();
    Mat3d(int n, int m, int k);
    inline T** operator[](const int i);//subscripting: pointer to row i
    inline const T* const * operator[](const int i) const;
    inline int dim1() const;
    inline int dim2() const;
    inline int dim3() const;
    ~Mat3d();
};

template <class T>
Mat3d<T>::Mat3d(): nn(0), mm(0), kk(0), v(0) {}

template <class T>
Mat3d<T>::Mat3d(int n, int m, int k) : nn(n), mm(m), kk(k), v(new T**[n])
{
    int i,j;
    v[0] = new T*[n*m];
    v[0][0] = new T[n*m*k];
    for(j=1; j<m; j++)
        v[0][j] = v[0][j-1] + k;
    for(i=1; i<n; i++) {
        v[i] = v[i-1] + m;
        v[i][0] = v[i-1][0] + m*k;
        for(j=1; j<m; j++)
            v[i][j] = v[i][j-1] + k;
    }
}

template <class T>
inline T** Mat3d<T>::operator[](const int i) //subscripting: pointer to row i
{
    return v[i];
}

template <class T>
inline const T* const * Mat3d<T>::operator[](const int i) const
{
    return v[i];
}

template <class T>
inline int Mat3d<T>::dim1() const
{
    return nn;
}

template <class T>
inline int Mat3d<T>::dim2() const
{
    return mm;
}

template <class T>
inline int Mat3d<T>::dim3() const
{
    return kk;
}

template <class T>
Mat3d<T>::~Mat3d()
{
    if (v != 0) {
        delete[] (v[0][0]);
        delete[] (v[0]);
        delete[] (v);
    }
}

typedef const Mat3d<double> Mat3D_I_DP;
typedef Mat3d<double> Mat3D_DP, Mat3D_O_DP, Mat3D_IO_DP;

#endif // __Mat3D_h__
